#ifndef EDDA_READER
#define EDDA_READER

#include <string>
#include "edda_export.h"
#include "dataset/dataset.h"

namespace edda{
	std::shared_ptr<Dataset<Real> > EDDA_EXPORT loadEddaScalarDataset_noneVTK(const std::string &edda_file);
	std::shared_ptr<Dataset<VECTOR3> > EDDA_EXPORT loadEddaVector3Dataset_noneVTK(const std::string &edda_file);
}

#endif
