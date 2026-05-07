// Minimal HDF5 contract checker for GR1D restart files consumed by the
// GR1DRestartHDF5 VisIt reader.

#include <hdf5.h>

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::vector<hsize_t> dataset_dims(hid_t file, const std::string& path) {
  hid_t ds = H5Dopen(file, path.c_str(), H5P_DEFAULT);
  if (ds < 0) throw std::runtime_error("missing dataset: " + path);

  hid_t space = H5Dget_space(ds);
  if (space < 0) {
    H5Dclose(ds);
    throw std::runtime_error("missing dataspace: " + path);
  }

  const int rank = H5Sget_simple_extent_ndims(space);
  if (rank <= 0) {
    H5Sclose(space);
    H5Dclose(ds);
    throw std::runtime_error("invalid rank for: " + path);
  }

  std::vector<hsize_t> dims(static_cast<std::size_t>(rank), 0);
  H5Sget_simple_extent_dims(space, dims.data(), nullptr);
  H5Sclose(space);
  H5Dclose(ds);
  return dims;
}

int read_i32(hid_t file, const std::string& path) {
  hid_t ds = H5Dopen(file, path.c_str(), H5P_DEFAULT);
  if (ds < 0) throw std::runtime_error("missing dataset: " + path);
  int value = 0;
  if (H5Dread(ds, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value) < 0) {
    H5Dclose(ds);
    throw std::runtime_error("failed to read: " + path);
  }
  H5Dclose(ds);
  return value;
}

double read_f64(hid_t file, const std::string& path) {
  hid_t ds = H5Dopen(file, path.c_str(), H5P_DEFAULT);
  if (ds < 0) throw std::runtime_error("missing dataset: " + path);
  double value = 0.0;
  if (H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value) < 0) {
    H5Dclose(ds);
    throw std::runtime_error("failed to read: " + path);
  }
  H5Dclose(ds);
  return value;
}

std::vector<double> read_f64_1d(hid_t file, const std::string& path) {
  const std::vector<hsize_t> dims = dataset_dims(file, path);
  if (dims.size() != 1) throw std::runtime_error(path + " must be rank 1");

  std::vector<double> values(static_cast<std::size_t>(dims[0]), 0.0);
  hid_t ds = H5Dopen(file, path.c_str(), H5P_DEFAULT);
  if (H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()) < 0) {
    H5Dclose(ds);
    throw std::runtime_error("failed to read: " + path);
  }
  H5Dclose(ds);
  return values;
}

void require_shape(hid_t file, const std::string& path, const std::vector<hsize_t>& expected) {
  const std::vector<hsize_t> actual = dataset_dims(file, path);
  if (actual != expected) throw std::runtime_error("unexpected shape for: " + path);
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 2) {
      std::cerr << "usage: gr1d_visit_reader_stub <restart.h5>\n";
      return 2;
    }

    H5Eset_auto(H5E_DEFAULT, nullptr, nullptr);
    hid_t file = H5Fopen(argv[1], H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) throw std::runtime_error("failed to open file");

    const int n1 = read_i32(file, "/n1");
    const int nt = read_i32(file, "/nt");
    const double time = read_f64(file, "/time");
    const std::vector<double> x1 = read_f64_1d(file, "/x1");
    if (static_cast<int>(x1.size()) != n1) throw std::runtime_error("x1 length != n1");

    int ghosts = 0;
    while (ghosts < n1 && x1[static_cast<std::size_t>(ghosts)] < 0.0) ++ghosts;
    const int interior = n1 - 2 * ghosts;
    if (interior <= 0) throw std::runtime_error("invalid inferred ghost count");

    require_shape(file, "/rho", {static_cast<hsize_t>(n1)});
    require_shape(file, "/press", {static_cast<hsize_t>(n1)});
    require_shape(file, "/temperature", {static_cast<hsize_t>(n1)});

    if (H5Lexists(file, "/q_M1", H5P_DEFAULT) > 0) {
      const int ng = read_i32(file, "/number_groups");
      const int ns = read_i32(file, "/number_species");
      require_shape(file, "/q_M1", {3, static_cast<hsize_t>(ng), static_cast<hsize_t>(ns),
                                    static_cast<hsize_t>(n1)});
    }

    H5Fclose(file);
    std::cout << "GR1D restart OK: n1=" << n1 << " ghosts=" << ghosts
              << " interior=" << interior << " nt=" << nt << " time=" << time << "\n";
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "gr1d_visit_reader_stub: " << ex.what() << "\n";
    return 3;
  }
}
