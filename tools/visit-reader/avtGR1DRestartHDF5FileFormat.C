// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

#include "avtGR1DRestartHDF5FileFormat.h"

#include <avtDatabaseMetaData.h>
#include <DebugStream.h>
#include <InvalidFilesException.h>
#include <InvalidVariableException.h>

#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkRectilinearGrid.h>

#include <algorithm>
#include <cctype>
#include <cstring>
#include <dirent.h>
#include <iomanip>
#include <regex>
#include <sstream>

using std::string;
using std::vector;

namespace
{
bool
read_all_double(hid_t dataset, vector<double> &values)
{
    hid_t dataspace = H5Dget_space(dataset);
    if (dataspace < 0)
        return false;

    const int rank = H5Sget_simple_extent_ndims(dataspace);
    if (rank <= 0)
    {
        H5Sclose(dataspace);
        return false;
    }

    vector<hsize_t> dims(static_cast<size_t>(rank), 0);
    H5Sget_simple_extent_dims(dataspace, dims.data(), NULL);

    size_t nelem = 1;
    for (size_t i = 0; i < dims.size(); ++i)
        nelem *= static_cast<size_t>(dims[i]);

    values.assign(nelem, 0.0);
    const herr_t status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                                  H5P_DEFAULT, values.data());
    H5Sclose(dataspace);
    return status >= 0;
}
}

avtGR1DRestartHDF5FileFormat::avtGR1DRestartHDF5FileFormat(const char *fname)
    : avtMTMDFileFormat(fname), initialized(false), n1(0), ghostCount(0),
      interiorStart(0), interiorCount(0)
{
}

avtGR1DRestartHDF5FileFormat::~avtGR1DRestartHDF5FileFormat()
{
    FreeUpResources();
}

void
avtGR1DRestartHDF5FileFormat::Initialize()
{
    if (initialized)
        return;

    H5Eset_auto(H5E_DEFAULT, NULL, NULL);
    DiscoverFiles();
    ReadFirstFileMetadata();
    initialized = true;
}

void
avtGR1DRestartHDF5FileFormat::DiscoverFiles()
{
    files.clear();

    const string selected(filename);
    const string directory = DirectoryName(selected);
    const string base = BaseName(selected);

    if (LooksLikeRestartName(base))
    {
        DIR *dir = opendir(directory.c_str());
        if (dir != NULL)
        {
            struct dirent *entry = NULL;
            while ((entry = readdir(dir)) != NULL)
            {
                const string candidateBase(entry->d_name);
                if (!LooksLikeRestartName(candidateBase))
                    continue;

                const string candidatePath = directory + "/" + candidateBase;
                int cycle = 0;
                double time = 0.0;
                if (IsRestartFile(candidatePath, &cycle, &time))
                {
                    FileState state;
                    state.path = candidatePath;
                    state.cycle = cycle;
                    state.time = time;
                    files.push_back(state);
                }
            }
            closedir(dir);
        }
    }

    if (files.empty())
    {
        int cycle = 0;
        double time = 0.0;
        if (IsRestartFile(selected, &cycle, &time))
        {
            FileState state;
            state.path = selected;
            state.cycle = cycle;
            state.time = time;
            files.push_back(state);
        }
    }

    if (files.empty())
        EXCEPTION1(InvalidFilesException, filename);

    std::sort(files.begin(), files.end(), [](const FileState &a, const FileState &b) {
        if (a.cycle != b.cycle)
            return a.cycle < b.cycle;
        return a.path < b.path;
    });
}

void
avtGR1DRestartHDF5FileFormat::ReadFirstFileMetadata()
{
    hid_t file = H5Fopen(files[0].path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0)
        EXCEPTION1(InvalidFilesException, files[0].path.c_str());

    if (!ReadScalarInt(file, "/n1", n1))
    {
        vector<hsize_t> dims;
        if (!ReadDatasetShape(file, "/x1", dims) || dims.size() != 1)
        {
            H5Fclose(file);
            EXCEPTION1(InvalidFilesException, files[0].path.c_str());
        }
        n1 = static_cast<int>(dims[0]);
    }

    vector<double> x = ReadX(file);
    if (x.size() != static_cast<size_t>(n1))
    {
        H5Fclose(file);
        EXCEPTION1(InvalidFilesException, files[0].path.c_str());
    }

    ghostCount = 0;
    while (ghostCount < n1 && x[static_cast<size_t>(ghostCount)] < 0.0)
        ++ghostCount;

    interiorStart = ghostCount;
    interiorCount = n1 - 2 * ghostCount;
    if (interiorCount <= 0)
    {
        interiorStart = 0;
        interiorCount = n1;
        ghostCount = 0;
    }

    interiorX = CropInterior(x);
    DiscoverVariables(file);
    H5Fclose(file);
}

bool
avtGR1DRestartHDF5FileFormat::IsRestartFile(const string &path, int *cycle,
                                            double *time) const
{
    hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0)
        return false;

    int localCycle = 0;
    double localTime = 0.0;
    vector<hsize_t> xdims;
    const bool ok = DatasetExists(file, "/x1") &&
                    DatasetExists(file, "/n1") &&
                    ReadScalarInt(file, "/nt", localCycle) &&
                    ReadScalarDouble(file, "/time", localTime) &&
                    ReadDatasetShape(file, "/x1", xdims) &&
                    xdims.size() == 1 && xdims[0] > 0;

    H5Fclose(file);
    if (!ok)
        return false;

    if (cycle != NULL)
        *cycle = localCycle;
    if (time != NULL)
        *time = localTime;
    return true;
}

void
avtGR1DRestartHDF5FileFormat::DiscoverVariables(hid_t file)
{
    variables.clear();
    variableMap.clear();

    hid_t root = H5Gopen(file, "/", H5P_DEFAULT);
    if (root < 0)
        return;

    H5G_info_t groupInfo;
    H5Gget_info(root, &groupInfo);
    for (hsize_t i = 0; i < groupInfo.nlinks; ++i)
    {
        char nameBuffer[512];
        const ssize_t n = H5Lget_name_by_idx(root, ".", H5_INDEX_NAME, H5_ITER_NATIVE,
                                             i, nameBuffer, sizeof(nameBuffer), H5P_DEFAULT);
        if (n <= 0)
            continue;

        const string name(nameBuffer);
        if (name == "x1")
            continue;

        const string path = "/" + name;
        vector<hsize_t> dims;
        if (!ReadDatasetShape(file, path, dims))
            continue;

        hid_t dataset = H5Dopen(file, path.c_str(), H5P_DEFAULT);
        if (dataset < 0)
            continue;
        const bool numeric = IsNumericDataset(dataset);
        H5Dclose(dataset);
        if (!numeric)
            continue;

        if (dims.size() == 1 && dims[0] == static_cast<hsize_t>(n1))
        {
            VarInfo info;
            info.name = name;
            info.dataset = path;
            info.kind = VarRank1;
            info.component = info.moment = info.group = info.species = 0;
            variables.push_back(info);
            variableMap[info.name] = info;
        }
        else if (name == "M1_matter_source" && dims.size() == 2 &&
                 dims[1] == static_cast<hsize_t>(n1))
        {
            for (int c = 0; c < static_cast<int>(dims[0]); ++c)
            {
                VarInfo info;
                std::ostringstream oss;
                oss << "M1_matter_source_c" << c + 1;
                info.name = oss.str();
                info.dataset = path;
                info.kind = VarM1MatterSource;
                info.component = c;
                info.moment = info.group = info.species = 0;
                variables.push_back(info);
                variableMap[info.name] = info;
            }
        }
        else if ((name == "q_M1" || name == "q_M1_fluid" || name == "eas") &&
                 dims.size() == 4 && dims[3] == static_cast<hsize_t>(n1))
        {
            const bool isEas = (name == "eas");
            for (int moment = 0; moment < static_cast<int>(dims[0]); ++moment)
            {
                for (int group = 0; group < static_cast<int>(dims[1]); ++group)
                {
                    for (int species = 0; species < static_cast<int>(dims[2]); ++species)
                    {
                        VarInfo info;
                        info.name = IndexedName(name, species + 1, group + 1,
                                                moment + 1, isEas ? "k" : "m");
                        info.dataset = path;
                        info.kind = isEas ? VarEAS :
                                    ((name == "q_M1") ? VarQ : VarQFluid);
                        info.component = 0;
                        info.moment = moment;
                        info.group = group;
                        info.species = species;
                        variables.push_back(info);
                        variableMap[info.name] = info;
                    }
                }
            }
        }
    }

    H5Gclose(root);
}

int
avtGR1DRestartHDF5FileFormat::GetNTimesteps(void)
{
    if (!initialized)
        Initialize();
    return static_cast<int>(files.size());
}

void
avtGR1DRestartHDF5FileFormat::GetCycles(vector<int> &cycles)
{
    if (!initialized)
        Initialize();
    cycles.clear();
    for (size_t i = 0; i < files.size(); ++i)
        cycles.push_back(files[i].cycle);
}

void
avtGR1DRestartHDF5FileFormat::GetTimes(vector<double> &times)
{
    if (!initialized)
        Initialize();
    times.clear();
    for (size_t i = 0; i < files.size(); ++i)
        times.push_back(files[i].time);
}

void
avtGR1DRestartHDF5FileFormat::FreeUpResources(void)
{
}

void
avtGR1DRestartHDF5FileFormat::PopulateDatabaseMetaData(avtDatabaseMetaData *md,
                                                       int timeState)
{
    (void)timeState;
    if (!initialized)
        Initialize();

    for (size_t i = 0; i < variables.size(); ++i)
    {
        avtCurveMetaData *cmd = new avtCurveMetaData;
        cmd->name = variables[i].name;
        cmd->xLabel = "radius";
        cmd->yLabel = variables[i].name;
        cmd->hasUnits = false;
        md->Add(cmd);
    }
}

vtkDataSet *
avtGR1DRestartHDF5FileFormat::GetMesh(int timeState, int domain,
                                      const char *meshname)
{
    (void)domain;
    if (!initialized)
        Initialize();

    const string varName(meshname);
    std::map<string, VarInfo>::const_iterator it = variableMap.find(varName);
    if (it == variableMap.end())
        EXCEPTION1(InvalidVariableException, meshname);

    hid_t file = OpenFile(timeState);
    if (file < 0)
        EXCEPTION1(InvalidFilesException, files[0].path.c_str());

    vector<double> data = ReadVariable(file, it->second);
    CloseFile(file);
    return BuildCurve(varName, interiorX, data);
}

vtkDataArray *
avtGR1DRestartHDF5FileFormat::GetVar(int timeState, int domain,
                                     const char *varname)
{
    (void)domain;
    if (!initialized)
        Initialize();

    std::map<string, VarInfo>::const_iterator it = variableMap.find(string(varname));
    if (it == variableMap.end())
        EXCEPTION1(InvalidVariableException, varname);

    hid_t file = OpenFile(timeState);
    if (file < 0)
        EXCEPTION1(InvalidFilesException, files[0].path.c_str());

    vector<double> data = ReadVariable(file, it->second);
    CloseFile(file);

    vtkFloatArray *arr = vtkFloatArray::New();
    arr->SetNumberOfTuples(static_cast<vtkIdType>(data.size()));
    arr->SetName(varname);
    for (size_t i = 0; i < data.size(); ++i)
        arr->SetValue(static_cast<vtkIdType>(i), static_cast<float>(data[i]));
    return arr;
}

vtkDataArray *
avtGR1DRestartHDF5FileFormat::GetVectorVar(int, int, const char *)
{
    return NULL;
}

hid_t
avtGR1DRestartHDF5FileFormat::OpenFile(int timeState)
{
    if (timeState < 0 || timeState >= static_cast<int>(files.size()))
        return -1;
    return H5Fopen(files[static_cast<size_t>(timeState)].path.c_str(),
                   H5F_ACC_RDONLY, H5P_DEFAULT);
}

void
avtGR1DRestartHDF5FileFormat::CloseFile(hid_t file)
{
    if (file >= 0)
        H5Fclose(file);
}

bool
avtGR1DRestartHDF5FileFormat::DatasetExists(hid_t file, const string &path) const
{
    hid_t dataset = H5Dopen(file, path.c_str(), H5P_DEFAULT);
    if (dataset < 0)
        return false;
    H5Dclose(dataset);
    return true;
}

bool
avtGR1DRestartHDF5FileFormat::ReadDatasetShape(hid_t file, const string &path,
                                               vector<hsize_t> &dims) const
{
    dims.clear();
    hid_t dataset = H5Dopen(file, path.c_str(), H5P_DEFAULT);
    if (dataset < 0)
        return false;

    hid_t dataspace = H5Dget_space(dataset);
    if (dataspace < 0)
    {
        H5Dclose(dataset);
        return false;
    }

    const int rank = H5Sget_simple_extent_ndims(dataspace);
    if (rank > 0)
    {
        dims.resize(static_cast<size_t>(rank), 0);
        H5Sget_simple_extent_dims(dataspace, dims.data(), NULL);
    }

    H5Sclose(dataspace);
    H5Dclose(dataset);
    return rank > 0;
}

bool
avtGR1DRestartHDF5FileFormat::ReadScalarInt(hid_t file, const string &path,
                                            int &value) const
{
    hid_t dataset = H5Dopen(file, path.c_str(), H5P_DEFAULT);
    if (dataset < 0)
        return false;
    const herr_t status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                                  H5P_DEFAULT, &value);
    H5Dclose(dataset);
    return status >= 0;
}

bool
avtGR1DRestartHDF5FileFormat::ReadScalarDouble(hid_t file, const string &path,
                                               double &value) const
{
    hid_t dataset = H5Dopen(file, path.c_str(), H5P_DEFAULT);
    if (dataset < 0)
        return false;
    const herr_t status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                                  H5P_DEFAULT, &value);
    H5Dclose(dataset);
    return status >= 0;
}

bool
avtGR1DRestartHDF5FileFormat::IsNumericDataset(hid_t dataset) const
{
    hid_t type = H5Dget_type(dataset);
    if (type < 0)
        return false;
    const H5T_class_t klass = H5Tget_class(type);
    H5Tclose(type);
    return klass == H5T_INTEGER || klass == H5T_FLOAT;
}

vector<double>
avtGR1DRestartHDF5FileFormat::ReadX(hid_t file)
{
    return ReadRank1(file, "/x1");
}

vector<double>
avtGR1DRestartHDF5FileFormat::ReadVariable(hid_t file, const VarInfo &info)
{
    switch (info.kind)
    {
      case VarRank1:
        return CropInterior(ReadRank1(file, info.dataset));
      case VarM1MatterSource:
        return ReadM1MatterSource(file, info);
      case VarQ:
      case VarQFluid:
      case VarEAS:
        return ReadRank4Slice(file, info);
    }
    return vector<double>();
}

vector<double>
avtGR1DRestartHDF5FileFormat::ReadRank1(hid_t file, const string &path)
{
    vector<double> values;
    hid_t dataset = H5Dopen(file, path.c_str(), H5P_DEFAULT);
    if (dataset < 0)
        return values;
    read_all_double(dataset, values);
    H5Dclose(dataset);
    return values;
}

vector<double>
avtGR1DRestartHDF5FileFormat::ReadM1MatterSource(hid_t file, const VarInfo &info)
{
    vector<double> values(static_cast<size_t>(interiorCount), 0.0);
    hid_t dataset = H5Dopen(file, info.dataset.c_str(), H5P_DEFAULT);
    if (dataset < 0)
        return values;

    hid_t dataspace = H5Dget_space(dataset);
    hsize_t start[2] = {static_cast<hsize_t>(info.component),
                        static_cast<hsize_t>(interiorStart)};
    hsize_t count[2] = {1, static_cast<hsize_t>(interiorCount)};
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, NULL, count, NULL);

    hsize_t memDims[1] = {static_cast<hsize_t>(interiorCount)};
    hid_t memspace = H5Screate_simple(1, memDims, NULL);
    H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT,
            values.data());

    H5Sclose(memspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);
    return values;
}

vector<double>
avtGR1DRestartHDF5FileFormat::ReadRank4Slice(hid_t file, const VarInfo &info)
{
    vector<double> values(static_cast<size_t>(interiorCount), 0.0);
    hid_t dataset = H5Dopen(file, info.dataset.c_str(), H5P_DEFAULT);
    if (dataset < 0)
        return values;

    hid_t dataspace = H5Dget_space(dataset);
    hsize_t start[4] = {static_cast<hsize_t>(info.moment),
                        static_cast<hsize_t>(info.group),
                        static_cast<hsize_t>(info.species),
                        static_cast<hsize_t>(interiorStart)};
    hsize_t count[4] = {1, 1, 1, static_cast<hsize_t>(interiorCount)};
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, NULL, count, NULL);

    hsize_t memDims[1] = {static_cast<hsize_t>(interiorCount)};
    hid_t memspace = H5Screate_simple(1, memDims, NULL);
    H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT,
            values.data());

    H5Sclose(memspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);
    return values;
}

vector<double>
avtGR1DRestartHDF5FileFormat::CropInterior(const vector<double> &values) const
{
    if (static_cast<int>(values.size()) < interiorStart + interiorCount)
        return vector<double>();
    return vector<double>(values.begin() + interiorStart,
                          values.begin() + interiorStart + interiorCount);
}

vtkDataSet *
avtGR1DRestartHDF5FileFormat::BuildCurve(const string &name,
                                         const vector<double> &x,
                                         const vector<double> &y)
{
    const size_t n = std::min(x.size(), y.size());
    vtkRectilinearGrid *grid = vtkRectilinearGrid::New();
    grid->SetDimensions(static_cast<int>(n), 1, 1);

    vtkFloatArray *xCoords = vtkFloatArray::New();
    xCoords->SetNumberOfTuples(static_cast<vtkIdType>(n));
    for (size_t i = 0; i < n; ++i)
        xCoords->SetValue(static_cast<vtkIdType>(i), static_cast<float>(x[i]));
    grid->SetXCoordinates(xCoords);
    xCoords->Delete();

    vtkFloatArray *yCoords = vtkFloatArray::New();
    yCoords->SetNumberOfTuples(1);
    yCoords->SetValue(0, 0.0f);
    grid->SetYCoordinates(yCoords);
    yCoords->Delete();

    vtkFloatArray *zCoords = vtkFloatArray::New();
    zCoords->SetNumberOfTuples(1);
    zCoords->SetValue(0, 0.0f);
    grid->SetZCoordinates(zCoords);
    zCoords->Delete();

    vtkFloatArray *values = vtkFloatArray::New();
    values->SetName(name.c_str());
    values->SetNumberOfTuples(static_cast<vtkIdType>(n));
    for (size_t i = 0; i < n; ++i)
        values->SetValue(static_cast<vtkIdType>(i), static_cast<float>(y[i]));
    grid->GetPointData()->SetScalars(values);
    values->Delete();

    return grid;
}

string
avtGR1DRestartHDF5FileFormat::DirectoryName(const string &path) const
{
    const size_t pos = path.rfind('/');
    if (pos == string::npos)
        return ".";
    if (pos == 0)
        return "/";
    return path.substr(0, pos);
}

string
avtGR1DRestartHDF5FileFormat::BaseName(const string &path) const
{
    const size_t pos = path.rfind('/');
    if (pos == string::npos)
        return path;
    return path.substr(pos + 1);
}

bool
avtGR1DRestartHDF5FileFormat::LooksLikeRestartName(const string &base) const
{
    static const std::regex pattern("^restart(_nt_.*_time_.*|_.*)?\\.h5$");
    return std::regex_match(base, pattern);
}

string
avtGR1DRestartHDF5FileFormat::IndexedName(const string &base, int species,
                                          int group, int index,
                                          const string &indexPrefix) const
{
    std::ostringstream oss;
    oss << base << "_s" << species << "_g" << std::setw(2) << std::setfill('0')
        << group << "_" << indexPrefix << index;
    return oss.str();
}
