#pragma once
#include <imagina/evaluator.h>
#include <vector>

namespace PTWithCompression {
	using namespace Imagina;
	using real = SRReal;
	using complex = SRComplex;

	struct Waypoint {
		complex Z;
		ITUInt Iteration;
	};

	struct CompressedReference {
		std::vector<Waypoint> Waypoints;
		SRComplex C;
	};

	class ReferenceCompressor {
		CompressedReference &reference;
		complex z = 0.0, c = 0.0;
		ITUInt iteration = 0;
		real tolerance = 0x1p-32;

	public:
		ReferenceCompressor(CompressedReference &reference, SRComplex C, real tolerance = 0x1p-32)
			: reference(reference), c(C), tolerance(tolerance) {
			reference.Waypoints.clear();
			reference.C = C;
		}

		void Add(complex Z) {
			z = z * z + c;
			iteration++;
			if (chebyshev_norm(z - Z) > chebyshev_norm(Z) * tolerance) {
				reference.Waypoints.push_back({ Z, iteration });
				z = Z;
			}
		}

		void Finalize(complex Z) {
			iteration++;
			reference.Waypoints.push_back({ Z, iteration });
		}
	};

	class ReferenceDecompressor {
		using Iterator = std::vector<Waypoint>::const_iterator;

		const CompressedReference &reference;
		complex Z = 0.0;
		ITUInt iteration = 0;
		Iterator nextWaypoint;

	public:
		ReferenceDecompressor(const CompressedReference &reference) : reference(reference) {
			nextWaypoint = reference.Waypoints.begin();
		}

		complex Get() { return Z; }

		complex Next() {
			iteration++;
			if (iteration == nextWaypoint->Iteration) {
				Z = nextWaypoint->Z;
				nextWaypoint++;
			} else {
				Z = Z * Z + reference.C;
			}
			return Z;
		}

		bool End() {
			return nextWaypoint == reference.Waypoints.end();
		}

		complex Reset() {
			Z = 0.0;
			iteration = 0;
			nextWaypoint = reference.Waypoints.begin();
			return Z;
		}
	};

	class PTWithCompressionEvaluator {
		struct Output {
			double Value;
		};

		StandardEvaluationParameters parameters;

		CompressedReference reference;

	public:
		const PixelDataInfo *GetOutputInfo();

		void Prepare(const HPReal &x, const HPReal &y, HRReal radius, const StandardEvaluationParameters &parameters);
		void Evaluate(IRasterizingInterface rasterizingInterface);
	};
}

IMPLEMENT_INTERFACE(PTWithCompression::PTWithCompressionEvaluator, Imagina::IEvaluator);
