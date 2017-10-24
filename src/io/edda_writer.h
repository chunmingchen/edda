#ifndef EDDA_WRITER
#define EDDA_WRITER

#include <string>
#include "edda_export.h"
#include "dataset/dataset.h"

namespace edda{
	void EDDA_EXPORT writeEddaDataset(std::shared_ptr<Dataset<Real> >, const std::string &edda_file);
	void EDDA_EXPORT writeEddaDataset(std::shared_ptr<Dataset<VECTOR3> >, const std::string &edda_file);
}

#endif
