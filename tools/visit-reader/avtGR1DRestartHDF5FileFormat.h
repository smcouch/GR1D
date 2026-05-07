// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

#ifndef AVT_GR1D_RESTART_HDF5_FILE_FORMAT_H
#define AVT_GR1D_RESTART_HDF5_FILE_FORMAT_H

#include <avtMTMDFileFormat.h>

#include <hdf5.h>

#include <map>
#include <string>
#include <vector>

class avtGR1DRestartHDF5FileFormat : public avtMTMDFileFormat
{
  public:
                       avtGR1DRestartHDF5FileFormat(const char *);
    virtual           ~avtGR1DRestartHDF5FileFormat();

    virtual void       GetCycles(std::vector<int> &);
    virtual void       GetTimes(std::vector<double> &);
    virtual int        GetNTimesteps(void);

    virtual const char *GetType(void) { return "GR1DRestartHDF5"; }
    virtual void       FreeUpResources(void);

    virtual vtkDataSet   *GetMesh(int, int, const char *);
    virtual vtkDataArray *GetVar(int, int, const char *);
    virtual vtkDataArray *GetVectorVar(int, int, const char *);

  protected:
    virtual void       PopulateDatabaseMetaData(avtDatabaseMetaData *, int);

  private:
    enum VarKind
    {
        VarRank1,
        VarM1MatterSource,
        VarQ,
        VarQFluid,
        VarEAS
    };

    struct VarInfo
    {
        std::string name;
        std::string dataset;
        VarKind kind;
        int component;
        int moment;
        int group;
        int species;
    };

    struct FileState
    {
        std::string path;
        int cycle;
        double time;
    };

    void                   Initialize();
    void                   DiscoverFiles();
    void                   ReadFirstFileMetadata();
    void                   DiscoverVariables(hid_t);

    bool                   IsRestartFile(const std::string &, int *, double *) const;
    bool                   DatasetExists(hid_t, const std::string &) const;
    bool                   ReadDatasetShape(hid_t, const std::string &, std::vector<hsize_t> &) const;
    bool                   ReadScalarInt(hid_t, const std::string &, int &) const;
    bool                   ReadScalarDouble(hid_t, const std::string &, double &) const;
    bool                   IsNumericDataset(hid_t) const;

    hid_t                  OpenFile(int);
    void                   CloseFile(hid_t);

    std::vector<double>    ReadX(hid_t);
    std::vector<double>    ReadVariable(hid_t, const VarInfo &);
    std::vector<double>    ReadRank1(hid_t, const std::string &);
    std::vector<double>    ReadM1MatterSource(hid_t, const VarInfo &);
    std::vector<double>    ReadRank4Slice(hid_t, const VarInfo &);
    std::vector<double>    CropInterior(const std::vector<double> &) const;
    vtkDataSet            *BuildCurve(const std::string &, const std::vector<double> &,
                                      const std::vector<double> &);

    std::string            DirectoryName(const std::string &) const;
    std::string            BaseName(const std::string &) const;
    bool                   LooksLikeRestartName(const std::string &) const;
    std::string            IndexedName(const std::string &, int, int, int,
                                       const std::string &) const;

    bool                   initialized;
    int                    n1;
    int                    ghostCount;
    int                    interiorStart;
    int                    interiorCount;

    std::vector<FileState> files;
    std::vector<double>    interiorX;
    std::vector<VarInfo>   variables;
    std::map<std::string, VarInfo> variableMap;
};

#endif
