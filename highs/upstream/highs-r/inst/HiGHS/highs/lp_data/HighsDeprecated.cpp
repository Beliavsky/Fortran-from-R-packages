/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsDeprecated.cpp
 * @brief
 */
#include "HConfig.h"
#include "Highs.h"

HighsStatus Highs::setLogCallback(void (*user_log_callback)(HighsLogType,
                                                            const char*, void*),
                                  void* user_log_callback_data) {
  
  options_.log_options.user_log_callback = user_log_callback;
  options_.log_options.user_log_callback_data = user_log_callback_data;
  return HighsStatus::kOk;
}

HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                       const bool value) {
  
  return setOptionValue(option, value);
}

HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                       const HighsInt value) {
  
  return setOptionValue(option, value);
}

HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                       const double value) {
  
  return setOptionValue(option, value);
}

HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                       const std::string& value) {
  
  return setOptionValue(option, value);
}

HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                       const char* value) {
  
  return setOptionValue(option, value);
}

HighsStatus Highs::readHighsOptions(const std::string& filename) {
  
  return readOptions(filename);
}

HighsStatus Highs::passHighsOptions(const HighsOptions& options) {
  
  return passOptions(options);
}

HighsStatus Highs::getHighsOptionValue(const std::string& option, bool& value) {
  
  return getOptionValue(option, value);
}

HighsStatus Highs::getHighsOptionValue(const std::string& option,
                                       HighsInt& value) {
  
  return getOptionValue(option, value);
}

HighsStatus Highs::getHighsOptionValue(const std::string& option,
                                       double& value) {
  
  return getOptionValue(option, value);
}

HighsStatus Highs::getHighsOptionValue(const std::string& option,
                                       std::string& value) {
  
  return getOptionValue(option, value);
}

HighsStatus Highs::getHighsOptionType(const std::string& option,
                                      HighsOptionType& type) {
  
  return getOptionType(option, type);
}

HighsStatus Highs::resetHighsOptions() {
  
  return resetOptions();
}

HighsStatus Highs::writeHighsOptions(
    const std::string& filename, const bool report_only_non_default_values) {
  
  return writeOptions(filename, report_only_non_default_values);
}

const HighsOptions& Highs::getHighsOptions() const {
  
  return getOptions();
}

HighsStatus Highs::setHighsLogfile(FILE* logfile) {
  
  options_.output_flag = false;
  return HighsStatus::kOk;
}

HighsStatus Highs::setHighsOutput(FILE* output) {
  
  options_.output_flag = false;
  return HighsStatus::kOk;
}

const HighsInfo& Highs::getHighsInfo() const {
  
  return getInfo();
}

HighsStatus Highs::getHighsInfoValue(const std::string& info, HighsInt& value) {
  
  return getInfoValue(info, value);
}

HighsStatus Highs::getHighsInfoValue(const std::string& info,
                                     double& value) const {
  
  return getInfoValue(info, value);
}

HighsStatus Highs::writeHighsInfo(const std::string& filename) {
  
  return writeInfo(filename);
}

double Highs::getHighsInfinity() {
  
  return getInfinity();
}

double Highs::getHighsRunTime() {
  
  return getRunTime();
}

#if 0
HighsStatus Highs::writeSolution(const std::string& filename,
                                const bool pretty) const {
  
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  FILE* file;
  HighsFileType file_type;
  call_status = openWriteFile(filename, "writeSolution", file, file_type);
  return_status =
      interpretCallStatus(call_status, return_status, "openWriteFile");
  if (return_status == HighsStatus::kError) return return_status;
  HighsInt style;
  if (pretty) {
    style = kSolutionStylePretty;
  } else {
    style = kSolutionStyleRaw;
  }
  writeSolutionFile(file, options_,
		    model_, basis_, solution_, info_, model_status_,
                    style);
  if (1) fclose(file);
  return HighsStatus::kOk;
}
#endif

const HighsModelStatus& Highs::getModelStatus(const bool) const {
  deprecationMessage("getModelStatus(const bool scaled_model)",
                     "getModelStatus()");
  return model_status_;
}
