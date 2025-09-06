#pragma once
#include <imagina/evaluator.h>
#include <vector>

namespace HarmonicLLA {
	using namespace Imagina;

	inline real_sr magnitude(const complex_sr &z) {
		//return std::abs(z);
		return chebyshev_norm(z);
	}

	struct LAStep {
		using real = real_sr;
		using complex = complex_sr;

		static constexpr real ValidRadiusScale = 0x1.0p-24;
		static constexpr real DipDetectionThreshold = 0x1.0p-10;

		complex Z, A, B;
		real ValidRadius, ValidRadiusC;
		size_t Length, NextStageLAIndex;

		LAStep() = default;
		LAStep(const LAStep &) = default;
		LAStep(LAStep &&) = default;
		explicit LAStep(size_t i, complex z) :
			Z(z), A(1.0), B(1.0), ValidRadius(1.0), ValidRadiusC(1.0), Length(1), NextStageLAIndex(i) {}

		LAStep(size_t i, complex z0, complex z1) :
			Z(z0), A(real(2.0) * z1), B(A + 1.0),
			ValidRadius(magnitude(z1) * ValidRadiusScale),
			ValidRadiusC(ValidRadius),
			Length(2), NextStageLAIndex(i) {}

		LAStep &operator=(const LAStep &) = default;
		LAStep &operator=(LAStep &&) = default;

		bool DetectDip(complex z) const {
			return magnitude(z) * ValidRadiusScale / magnitude(A) < ValidRadius * DipDetectionThreshold;
		}

		std::pair<LAStep, bool> Step(complex z) const {
			LAStep result;
			real radius = magnitude(z) * ValidRadiusScale;

			result.ValidRadius = std::min(ValidRadius, radius / magnitude(A));
			result.ValidRadiusC = std::min(ValidRadiusC, radius / magnitude(B));
			result.A = 2.0 * z * A;
			result.B = 2.0 * z * B + 1.0;
			result.Z = Z;

			result.Length = Length + 1;
			result.NextStageLAIndex = NextStageLAIndex;

			bool dipDetected = result.ValidRadius < ValidRadius * DipDetectionThreshold;

			return { result, dipDetected };
		}

		std::pair<LAStep, bool> Composite(const LAStep &step) const {
			LAStep result;

			real radius = magnitude(step.Z) * ValidRadiusScale;

			result.ValidRadius = std::min(ValidRadius, radius / magnitude(A));
			result.ValidRadiusC = std::min(ValidRadiusC, radius / magnitude(B));
			result.A = 2.0 * step.Z * A;
			result.B = 2.0 * step.Z * B;

			bool dipDetected = result.ValidRadius < ValidRadius * DipDetectionThreshold;

			result.ValidRadius = std::min(result.ValidRadius, step.ValidRadius / magnitude(result.A));
			result.ValidRadiusC = std::min(result.ValidRadiusC, step.ValidRadius / magnitude(result.B));
			result.A = result.A * step.A;
			result.B = result.B * step.A + step.B;
			result.Z = Z;

			result.Length = Length + step.Length;
			result.NextStageLAIndex = NextStageLAIndex;

			return { result, dipDetected };
		}
	};

	class HarmonicLLAEvaluator {
		using real = real_sr;
		using complex = complex_sr;
		struct Output {
			double Value;
		};

		struct LAStageInfo {
			uint64_t Begin;
			uint64_t End;
		};

		StandardEvaluationParameters parameters;

		uint64_t referenceLength;
		complex *reference = nullptr;

		std::vector<LAStageInfo> LAStages;
		std::vector<LAStep> LASteps;

		void ComputeOrbit(const real_hp &x, const real_hp &y, real_hr radius);
		bool CreateLAFromOrbit();
		bool CreateNewLAStage();

	public:
		const PixelDataInfo *GetOutputInfo();

		void Prepare(const real_hp &x, const real_hp &y, real_hr radius, const StandardEvaluationParameters &parameters);
		void Evaluate(IRasterizer rasterizer);
	};
}

IMPLEMENT_INTERFACE(HarmonicLLA::HarmonicLLAEvaluator, Imagina::IEvaluator);
