#ifndef REPORT_DRUG_HPP
#define REPORT_DRUG_HPP

#include <map>
#include <string>
#include <types/drug_block_input.hpp>
#include <types/cipa_features.hpp>
#include <types/parameter.hpp>
#include <vector>

using std::multimap;
using std::string;
using std::vector;

int generate_report_drug(const Parameter *p_param, Cipa_Features *p_features, Drug_Block_Input &hill);

#endif
