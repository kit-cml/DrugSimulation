#ifndef REPORT_DRUG_HPP
#define REPORT_DRUG_HPP

#include <types/drug_block_input.hpp>
#include <types/parameter.hpp>

int generate_report_drug(const Parameter *p_param);

int get_concentrations_size(const char* concs_str);

#endif
