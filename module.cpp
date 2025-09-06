#include <imagina/module.h>
#include "Perturbation.h"
#include "PTWithCompression.h"
#include "HarmonicLLA.h"
#include "HarmonicMLA.h"
#include "MipLA.h"

using namespace Imagina;

constexpr ComponentInfo Components[]{
	ComponentInfo::Evaluator<Perturbation::PerturbationEvaluator>("PerturbationEvaluator", "Perturbation"),
	ComponentInfo::Evaluator<PTWithCompression::PTWithCompressionEvaluator>("PTWithCompression", "Perturbation With Compression"),
	ComponentInfo::Evaluator<HarmonicLLA::HarmonicLLAEvaluator>("HarmonicLLAEvaluator", "Harmonic Lin LA"),
	ComponentInfo::Evaluator<HarmonicMLA::HarmonicMLAEvaluator>("HarmonicMLAEvaluator", "Harmonic Mag LA"),
	ComponentInfo::Evaluator<MipLA::MipLAEvaluator>("MipLAEvaluator", "Mip LA"),
};
constexpr ModuleInfo Module("Algorithms", "Algorithms", Components);

im_api const Imagina::ModuleInfo *ImGetModuleInfo() {
	return &Module;
}