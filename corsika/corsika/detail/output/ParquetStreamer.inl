/*
 * (c) Copyright 2021 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

namespace corsika {

  inline ParquetStreamer::ParquetStreamer()
      : isInit_(false) {}

  inline void ParquetStreamer::initStreamer(std::string const& filepath) {

    // open the file and connect it to our pointer
    PARQUET_ASSIGN_OR_THROW(outfile_, arrow::io::FileOutputStream::Open(filepath));

    // the default builder settings
    builder_.created_by("CORSIKA8");

    // add run and event tags to the file
    addField("shower", parquet::Repetition::REQUIRED, parquet::Type::INT32,
             parquet::ConvertedType::UINT_32);
  }

  template <typename... TArgs>
  inline void ParquetStreamer::addField(TArgs&&... args) {
    fields_.push_back(parquet::schema::PrimitiveNode::Make(args...));
  }

  inline void ParquetStreamer::enableCompression(int const level) {
    builder_.compression(parquet::Compression::LZ4);
    builder_.compression_level(level);
  }

  inline void ParquetStreamer::buildStreamer() {

    // build the top level schema
    schema_ = std::static_pointer_cast<parquet::schema::GroupNode>(
        parquet::schema::GroupNode::Make("schema", parquet::Repetition::REQUIRED,
                                         fields_));

    // and build the writer
    writer_ = std::make_shared<parquet::StreamWriter>(
        parquet::ParquetFileWriter::Open(outfile_, schema_, builder_.build()));

    //  only now this object is ready to stream
    isInit_ = true;
  }

  inline void ParquetStreamer::closeStreamer() {
    writer_.reset();
    [[maybe_unused]] auto status = outfile_->Close();
    isInit_ = false;
  }

  inline std::shared_ptr<parquet::StreamWriter> ParquetStreamer::getWriter() {
    if (!isInit()) {
      throw std::runtime_error(
          "ParquetStreamer not initialized. Either 1) add the "
          "corresponding module to "
          "the OutputManager, or 2) declare the module to write no output using "
          "NoOutput.");
    }
    return writer_;
  }

} // namespace corsika
