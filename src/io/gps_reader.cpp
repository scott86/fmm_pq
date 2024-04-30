/**
 *  Implementation of GPS reader classes
 *
 * @author: Can Yang
 * @version: 2017.11.11
 */
#include "io/gps_reader.hpp"
#include "util/debug.hpp"
#include "util/util.hpp"
#include "config/gps_config.hpp"
#include <iostream>
#include <stdexcept>
#include <string>

#include <boost/format.hpp>

#include "arrow/io/api.h"
#include <arrow/filesystem/api.h>
#include "parquet/arrow/reader.h"

using namespace FMM;
using namespace FMM::CORE;
using namespace FMM::IO;
using namespace FMM::UTIL;

std::vector<Trajectory> ITrajectoryReader::read_next_N_trajectories(int N) {
  std::vector<Trajectory> trajectories;
  int i = 0;
  while (i < N && has_next_trajectory()) {
    trajectories.push_back(read_next_trajectory());
    ++i;
  }
  return trajectories;
}

std::vector<Trajectory> ITrajectoryReader::read_all_trajectories() {
  std::vector<Trajectory> trajectories;
  int i = 0;
  while (has_next_trajectory()) {
    trajectories.push_back(read_next_trajectory());
    ++i;
  }
  return trajectories;
}

GDALTrajectoryReader::GDALTrajectoryReader(const std::string &filename,
                                           const std::string &id_name,
                                           const std::string &timestamp_name) {
  SPDLOG_INFO("Read trajectory from file {}",filename);
  OGRRegisterAll();
  poDS = (GDALDataset *) GDALOpenEx(filename.c_str(),
                                    GDAL_OF_VECTOR, NULL, NULL, NULL);
  if (poDS == NULL) {
    std::string message = "Open data source fail";
    SPDLOG_CRITICAL(message);
    throw std::runtime_error(message);
  }
  ogrlayer = poDS->GetLayer(0);
  _cursor = 0;
  // Get the number of features first
  OGRFeatureDefn *ogrFDefn = ogrlayer->GetLayerDefn();
  NUM_FEATURES = ogrlayer->GetFeatureCount();
  // This should be a local field rather than a new variable
  id_idx = ogrFDefn->GetFieldIndex(id_name.c_str());
  if (id_idx < 0) {
    std::string message = (boost::format("Id column %1% not found") % id_name).str();
    SPDLOG_CRITICAL(message);
    GDALClose(poDS);
    throw std::runtime_error(message);
  }
  timestamp_idx = ogrFDefn->GetFieldIndex(timestamp_name.c_str());
  if (timestamp_idx < 0) {
    SPDLOG_WARN("Timestamp column {} not found", timestamp_name);
  }
  if (wkbFlatten(ogrFDefn->GetGeomType()) != wkbLineString) {
    std::string message = (boost::format("Geometry type is %1%, which should be linestring") %
            OGRGeometryTypeToName(ogrFDefn->GetGeomType())).str();
    SPDLOG_CRITICAL(message);
    GDALClose(poDS);
    throw std::runtime_error(message);
  } else {
    SPDLOG_DEBUG("Geometry type is {}",
                 OGRGeometryTypeToName(ogrFDefn->GetGeomType()));
  }
  SPDLOG_INFO("Total number of trajectories {}", NUM_FEATURES);
  SPDLOG_INFO("Finish reading meta data");
}

bool GDALTrajectoryReader::has_next_trajectory() {
  return _cursor < NUM_FEATURES;
}

bool GDALTrajectoryReader::has_timestamp() {
  return timestamp_idx > 0;
}

Trajectory GDALTrajectoryReader::read_next_trajectory() {
  OGRFeature *ogrFeature = ogrlayer->GetNextFeature();
  int trid = ogrFeature->GetFieldAsInteger(id_idx);
  OGRGeometry *rawgeometry = ogrFeature->GetGeometryRef();
  FMM::CORE::LineString linestring =
    FMM::CORE::ogr2linestring((OGRLineString *) rawgeometry);
  OGRFeature::DestroyFeature(ogrFeature);
  ++_cursor;
  return Trajectory{trid, linestring};
}

int GDALTrajectoryReader::get_num_trajectories() {
  return NUM_FEATURES;
}

void GDALTrajectoryReader::close() {
  GDALClose(poDS);
}

CSVTrajectoryReader::CSVTrajectoryReader(const std::string &e_filename,
                                         const std::string &id_name,
                                         const std::string &geom_name,
                                         const std::string &timestamp_name) :
  ifs(e_filename) {
  std::string line;
  std::getline(ifs, line);
  std::stringstream check1(line);
  std::string intermediate;
  // Tokenizing w.r.t. space ' '
  int i = 0;
  while (safe_get_line(check1, intermediate, delim)) {
    if (intermediate == id_name) {
      id_idx = i;
    }
    if (intermediate == geom_name) {
      geom_idx = i;
    }
    if (intermediate == timestamp_name) {
      timestamp_idx = i;
    }
    ++i;
  }
  if (id_idx < 0 || geom_idx < 0) {
    std::string message = (boost::format("Id %1% or Geometry column %2% not found") % id_name % geom_name).str();
    SPDLOG_CRITICAL(message);
    throw std::runtime_error(message);
  }
  if (timestamp_idx < 0) {
    SPDLOG_WARN("Timestamp column {} not found", timestamp_name);
  }
  SPDLOG_INFO("Id index {} Geometry index {} Timstamp index {}",
              id_idx, geom_idx, timestamp_idx);
}

std::vector<double> CSVTrajectoryReader::string2time(
  const std::string &str) {
  std::vector<double> values;
  std::stringstream ss(str);
  double v;
  while (ss >> v) {
    values.push_back(v);
    if (ss.peek() == ',')
      ss.ignore();
  }
  return values;
}

bool CSVTrajectoryReader::has_timestamp() {
  return timestamp_idx > 0;
}

Trajectory CSVTrajectoryReader::read_next_trajectory() {
  // Read the geom idx column into a trajectory
  std::string line;
  std::getline(ifs, line);
  std::stringstream ss(line);
  int trid = 0;
  int index = 0;
  std::string intermediate;
  FMM::CORE::LineString geom;
  std::vector<double> timestamps;
  while (std::getline(ss, intermediate, delim)) {
    if (index == id_idx) {
      trid = std::stoi(intermediate);
    }
    if (index == geom_idx) {
      // intermediate
      boost::geometry::read_wkt(intermediate, geom.get_geometry());
    }
    if (index == timestamp_idx) {
      // intermediate
      timestamps = string2time(intermediate);
    }
    ++index;
  }
  return Trajectory{trid, geom, timestamps};
}

bool CSVTrajectoryReader::has_next_trajectory() {
  return ifs.peek() != EOF;
}

void CSVTrajectoryReader::reset_cursor() {
  ifs.clear();
  ifs.seekg(0, std::ios::beg);
  std::string line;
  std::getline(ifs, line);
}
void CSVTrajectoryReader::close() {
  ifs.close();
}

CSVPointReader::CSVPointReader(const std::string &e_filename,
                               const std::string &id_name,
                               const std::string &x_name,
                               const std::string &y_name,
                               const std::string &time_name) :
  ifs(e_filename) {
  std::string line;
  std::getline(ifs, line);
  std::stringstream check1(line);
  std::string intermediate;
  // Tokenizing w.r.t. space ' '
  int i = 0;
  while (safe_get_line(check1, intermediate, delim)) {
    if (intermediate == id_name) {
      id_idx = i;
    }
    if (intermediate == x_name) {
      x_idx = i;
    }
    if (intermediate == y_name) {
      y_idx = i;
    }
    if (intermediate == time_name) {
      timestamp_idx = i;
    }
    ++i;
  }
  if (id_idx < 0 || x_idx < 0 || y_idx < 0) {
    if (id_idx < 0) {
      std::string message = (boost::format("Id column %1% not found") % id_name).str();
      SPDLOG_CRITICAL(message);
      throw std::runtime_error(message);
    }
    if (x_idx < 0) {
      std::string message = (boost::format("X column name %1% not found") % x_name).str();
      SPDLOG_CRITICAL(message);
      throw std::runtime_error(message);
    }
    if (y_idx < 0) {
      std::string message = (boost::format("Y column name %1% not found") % y_name).str();
      SPDLOG_CRITICAL(message);
      throw std::runtime_error(message);
    }
  }
  if (timestamp_idx < 0) {
    SPDLOG_WARN("Time stamp {} not found, will be estimated ", time_name);
  }
  SPDLOG_INFO("Id index {} x index {} y index {} time index {}",
              id_idx, x_idx, y_idx, timestamp_idx);
}

Trajectory CSVPointReader::read_next_trajectory() {
  // Read the geom idx column into a trajectory
  std::string intermediate;
  FMM::CORE::LineString geom;
  std::vector<double> timestamps;
  bool on_same_trajectory = true;
  bool first_observation = true;
  int trid = -1;
  int prev_id = -1;
  double prev_timestamp = -1.0;
  std::string line;
  while (on_same_trajectory && has_next_trajectory()) {
    if (prev_line.empty()) {
      std::getline(ifs, line);
    } else {
      line = prev_line;
      prev_line.clear();
    }
    std::stringstream ss(line);
    int id = 0;
    double x = 0, y = 0;
    double timestamp = 0;
    int index = 0;
    while (std::getline(ss, intermediate, delim)) {
      if (index == id_idx) {
        id = std::stoi(intermediate);
      }
      if (index == x_idx) {
        x = std::stof(intermediate);
      }
      if (index == y_idx) {
        y = std::stof(intermediate);
      }
      if (index == timestamp_idx) {
        timestamp = std::stof(intermediate);
      }
      ++index;
    }
    if (prev_id == id || first_observation) {
      geom.add_point(x, y);
      if (has_timestamp())
        timestamps.push_back(timestamp);
    }
    if (prev_id != id && !first_observation) {
      on_same_trajectory = false;
      trid = prev_id;
    }
    first_observation = false;
    prev_id = id;
    if (!on_same_trajectory) {
      prev_line = line;
    }
  }
  if (!has_next_trajectory()) {
    trid = prev_id;
  }
  return Trajectory{trid, geom, timestamps};
}

bool CSVPointReader::has_next_trajectory() {
  return ifs.peek() != EOF;
}

void CSVPointReader::reset_cursor() {
  ifs.clear();
  ifs.seekg(0, std::ios::beg);
  std::string line;
  std::getline(ifs, line);
}

void CSVPointReader::close() {
  ifs.close();
}

bool CSVPointReader::has_timestamp() {
  return timestamp_idx > 0;
}


arrow::Result<std::unique_ptr<parquet::arrow::FileReader>> OpenReader(const std::string &e_filename) {
  arrow::fs::LocalFileSystem file_system;
  ARROW_ASSIGN_OR_RAISE(auto input, file_system.OpenInputFile(e_filename));

  std::cout << "open1" << std::endl;

  parquet::ArrowReaderProperties arrow_reader_properties =
      parquet::default_arrow_reader_properties();

  arrow_reader_properties.set_pre_buffer(true);
  arrow_reader_properties.set_use_threads(true);

  std::cout << "open2" << std::endl;

  parquet::ReaderProperties reader_properties =
      parquet::default_reader_properties();

  // Open Parquet file reader
  std::unique_ptr<parquet::arrow::FileReader> arrow_reader;
  auto reader_builder = parquet::arrow::FileReaderBuilder();
  reader_builder.properties(arrow_reader_properties);
  ARROW_RETURN_NOT_OK(reader_builder.Open(std::move(input), reader_properties));
  ARROW_RETURN_NOT_OK(reader_builder.Build(&arrow_reader));

  std::cout << "open3" << std::endl;

  return arrow_reader;
}


ParquetReader::ParquetReader(const std::string &e_filename) {
  std::cout << "debug1" << std::endl;

  // open file for reading
  init(e_filename);

  std::cout << "debug2" << std::endl;

  // initialize some stuff
  idx = 0;
  n = table->num_rows();
  curTripId = -1;

  std::cout << "debug3" << std::endl;

  std::cout << table->schema()->ToString() << std::endl;

  // get timestamp column deets
  timestampIdx      = 0;
  timestamp_col     = table->column(0);
  timestampChunks   = timestamp_col->num_chunks();
  curTimestampChunk = 0;
  timestamp_chunk   = std::dynamic_pointer_cast<arrow::TimestampArray>(timestamp_col->chunk(0));

  std::cout << "debug4" << std::endl;

  // get lon column deets
  lonIdx      = 0;
  lon_col     = table->column(2);
  lonChunks   = lon_col->num_chunks();
  curLonChunk = 0;
  lon_chunk   = std::dynamic_pointer_cast<arrow::DoubleArray>(lon_col->chunk(0));

  std::cout << "debug5" << std::endl;

  // get lat column deets
  latIdx      = 0;
  lat_col     = table->column(1);
  latChunks   = lat_col->num_chunks();
  curLatChunk = 0;
  lat_chunk   = std::dynamic_pointer_cast<arrow::DoubleArray>(lat_col->chunk(0));

  std::cout << "debug6" << std::endl;

  // get trip column deets
  tripIdx      = 0;
  trip_col     = table->column(3);
  tripChunks   = trip_col->num_chunks();
  curTripChunk = 0;
  trip_chunk   = std::dynamic_pointer_cast<arrow::Int32Array>(trip_col->chunk(0));

  std::cout << "debug7" << std::endl;
}

arrow::Status ParquetReader::init(const std::string &e_filename) {
  ARROW_ASSIGN_OR_RAISE(auto arrow_reader, OpenReader(e_filename));
  ARROW_RETURN_NOT_OK(arrow_reader->ReadTable(&table));
  return arrow::Status::OK();
}

// return current row and advance indicies to next spot
TripsRow ParquetReader::getCurrentRow() {

  // assemble current row
  //int agent     = agent_chunk->Value(agentIdx);
  double t      = static_cast<double>(timestamp_chunk->Value(timestampIdx)) / 1000.0; // millis -> seconds
  double lat    = lat_chunk->Value(latIdx);
  double lon    = lon_chunk->Value(lonIdx);
  int trip      = trip_chunk->Value(tripIdx);

  return TripsRow{t, lat, lon, trip};
}

// advance indicies
void ParquetReader::nextRow() {

  // advance master index
  idx = idx + 1;

  // don't bother continuing if we're at the end
  if(idx >= n) { return; }

  // advance each column
  nextTimestamp();
  nextLon();
  nextLat();
  nextTrip();
}

void ParquetReader::nextTimestamp() {
  timestampIdx = timestampIdx + 1;
  if(timestampIdx >= timestamp_chunk->length()) {
    // advance to next chunk
    timestampIdx = 0;
    curTimestampChunk = curTimestampChunk + 1;
    timestamp_chunk = std::dynamic_pointer_cast<arrow::TimestampArray>(timestamp_col->chunk(curTimestampChunk));
  }
}

void ParquetReader::nextLon() {
  lonIdx = lonIdx + 1;
  if(lonIdx >= lon_chunk->length()) {
    // advance to next chunk
    lonIdx = 0;
    curLonChunk = curLonChunk + 1;
    lon_chunk = std::dynamic_pointer_cast<arrow::DoubleArray>(lon_col->chunk(curLonChunk));
  }
}

void ParquetReader::nextLat() {
  latIdx = latIdx + 1;
  if(latIdx >= lat_chunk->length()) {
    // advance to next chunk
    latIdx = 0;
    curLatChunk = curLatChunk + 1;
    lat_chunk = std::dynamic_pointer_cast<arrow::DoubleArray>(lat_col->chunk(curLatChunk));
  }
}

void ParquetReader::nextTrip() {
  tripIdx = tripIdx + 1;
  if(tripIdx >= trip_chunk->length()) {
    // advance to next chunk
    tripIdx = 0;
    curTripChunk = curTripChunk + 1;
    trip_chunk = std::dynamic_pointer_cast<arrow::Int32Array>(trip_col->chunk(curTripChunk));
  }
}

bool ParquetReader::has_next_trajectory() {
  return idx < n;
}

Trajectory ParquetReader::read_next_trajectory() {
  //std::cout << "rn1" << std::endl;

  // initialize some vars
  FMM::CORE::LineString geom;
  std::vector<double> timestamps;

  //std::cout << "rn2" << std::endl;

  //std::cout << lon_chunk->length() << std::endl;

  // first leg
  //std::cout << tripIdx << std::endl;
  //std::cout << trip_chunk->length() << std::endl;
  curTripId = trip_chunk->Value(tripIdx);
  //std::cout << "rn2a" << std::endl;
  std::cout << timestamp_chunk->Value(timestampIdx) << std::endl;
  timestamps.push_back(timestamp_chunk->Value(timestampIdx));
  //std::cout << "rn2b" << std::endl;
  double x = lon_chunk->Value(lonIdx);
  //std::cout << "rn2c" << std::endl;
  double y = lat_chunk->Value(latIdx);
  //std::cout << "rn2d" << std::endl;
  geom.add_point(x, y);
  //std::cout << "rn2e" << std::endl;

  //std::cout << "rn3" << std::endl;

  int curTripLength = 1;
  
  // read rows until there is a change in trip ID
  while(has_next_trajectory()) {
    TripsRow data = getCurrentRow();
    nextRow();

    //std::cout << "rn4" << std::endl;

    // break if change in trip ID
    if(curTripId != data.trip) { break; }

    //std::cout << "rn5" << std::endl;

    // append current row to trip geom
    timestamps.push_back(data.timestamp);
    geom.add_point(data.lon, data.lat);
    curTripLength = curTripLength + 1;

    //std::cout << "rn6" << std::endl;
  }

  return Trajectory{curTripId, geom, timestamps};
}

// assume true for trips files
bool ParquetReader::has_timestamp() {
  return true;
}

// no resources to release?
void ParquetReader::close() {
}



GPSReader::GPSReader(const FMM::CONFIG::GPSConfig &config) {
  mode = config.get_gps_format();
  if (mode == 0) {
    SPDLOG_INFO("GPS data in trajectory shapefile format");
    reader = std::make_shared<GDALTrajectoryReader>
               (config.file, config.id,config.timestamp);
  } else if (mode == 1) {
    SPDLOG_INFO("GPS data in trajectory CSV format");
    reader = std::make_shared<CSVTrajectoryReader>
               (config.file, config.id, config.geom, config.timestamp);
  } else if (mode == 2) {
    SPDLOG_INFO("GPS data in point CSV format");
    reader = std::make_shared<CSVPointReader>
               (config.file, config.id, config.x, config.y, config.timestamp);
  } else if (mode == 3) {
    SPDLOG_INFO("GPS data in parquet format");
    reader = std::make_shared<ParquetReader>
               (config.file);

  } else {
    std::string message = "Unrecognized GPS format";
    SPDLOG_CRITICAL(message);
    throw std::runtime_error(message);
  }
};
