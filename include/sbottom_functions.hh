#ifndef SBOTTOM_FUNCTIONS_H
#define SBOTTOM_FUNCTIONS_H

class SUSYObjDef;
class SusyBuffer;

bool ChfCheck(const std::vector<size_t>& jet_indices, const SusyBuffer &buf, SUSYObjDef &def);

#endif