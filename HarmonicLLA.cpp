#include <imagina/pixel_management.h>
#include <imagina/output_info_helper.h>
#include "HarmonicLLA.h"

namespace HarmonicLLA {
	const Imagina::PixelDataInfo *HarmonicLLAEvaluator::GetOutputInfo() {
		IM_GET_OUTPUT_INFO_IMPL(Output, Value);
	}

	void HarmonicLLAEvaluator::Prepare(const real_hp &x, const real_hp &y, real_hr radius, const StandardEvaluationParameters &parameters) {
		this->parameters = parameters;
		delete[] reference;
		LAStages.clear();
		LASteps.clear();

		ComputeOrbit(x, y, radius);
		if (referenceLength < 8) return;

		if (!CreateLAFromOrbit()) return;
		while (CreateNewLAStage());
	}
	void HarmonicLLAEvaluator::ComputeOrbit(const real_hp &x, const real_hp &y, real_hr radius) {
		reference = new complex[parameters.Iterations + 1];

		complex_hp C = complex_hp(x, y);
		complex_hp Z = C;
		complex z = complex(Z), dzdc = 1.0;

		reference[0] = 0.0;
		reference[1] = z;

		size_t i = 1;
		while (i < parameters.Iterations) {
			dzdc = 2.0 * z * dzdc + 1.0;
			Z = Z * Z + C;
			z = complex(Z);
			i++;

			reference[i] = z;

			if (radius * chebyshev_norm(dzdc) * 2.0 > chebyshev_norm(z)) break;
			if (norm(z) > 16.0) break;
		};

		referenceLength = i;
	}

	bool HarmonicLLAEvaluator::CreateLAFromOrbit() {
		size_t Period = 0;

		LAStep step = LAStep(0, 0.0, reference[1]);

		size_t i;
		for (i = 2; i < referenceLength; i++) {
			auto [newStep, dipDetected] = step.Step(reference[i]);

			if (dipDetected) {
				Period = i;
				break;
			}
			step = newStep;
		}

		LASteps.push_back(step);

		if (!Period) {
			LAStages.push_back({ 0, 1 });
			LASteps.push_back(LAStep(0, reference[referenceLength]));

			return false;
		}

		if (i + 1 >= referenceLength) {
			step = LAStep(i, reference[i]);
			i++;
		} else {
			step = LAStep(i, reference[i], reference[i + 1]);
			i += 2;
		}
		for (; i < referenceLength; i++) {
			auto [newStep, dipDetected] = step.Step(reference[i]);

			if (dipDetected || step.Length >= Period) {
				LASteps.push_back(step);

				if (i + 1 >= referenceLength || newStep.DetectDip(reference[i + 1])) {
					step = LAStep(i, reference[i]);
				} else {
					step = LAStep(i, reference[i], reference[i + 1]);
					i++;
				}
			} else {
				step = newStep;
			}
		}

		LASteps.push_back(step);
		LAStages.push_back({ 0, LASteps.size() });
		LASteps.push_back(LAStep(0, reference[referenceLength]));

		return true;
	}

	bool HarmonicLLAEvaluator::CreateNewLAStage() {
		LAStageInfo prevStage = LAStages.back();
		size_t begin = LASteps.size();

		size_t Period = 0;
		size_t i = prevStage.Begin;

		LAStep step = LASteps[i].Composite(LASteps[i + 1]).first;
		step.NextStageLAIndex = i;
		i += 2;

		for (; i < prevStage.End; i++) {
			auto [newStep, dipDetected] = step.Composite(LASteps[i]);

			if (dipDetected) {
				Period = step.Length;

				LASteps.push_back(step);

				if (i + 1 >= prevStage.End || newStep.DetectDip(LASteps[i + 1].Z)) {
					step = LASteps[i];
					step.NextStageLAIndex = i;
					i++;
				} else {
					step = LASteps[i].Composite(LASteps[i + 1]).first;
					step.NextStageLAIndex = i;
					i += 2;
				}
				break;
			}
			step = newStep;
		}

		if (!Period) {
			LASteps.push_back(step);
			LAStages.push_back({ begin, LASteps.size() });
			LASteps.push_back(LASteps[prevStage.End]);

			return false;
		}

		for (; i < prevStage.End; i++) {
			auto [newStep, dipDetected] = step.Composite(LASteps[i]);

			if (dipDetected || step.Length >= Period) {
				LASteps.push_back(step);

				if (i + 1 >= prevStage.End || newStep.DetectDip(LASteps[i + 1].Z)) {
					step = LASteps[i];
					step.NextStageLAIndex = i;
				} else {
					step = LASteps[i].Composite(LASteps[i + 1]).first;
					step.NextStageLAIndex = i;
					i++;
				}
			} else {
				step = newStep;
			}
		}

		LASteps.push_back(step);
		LAStages.push_back({ begin, LASteps.size() });
		LASteps.push_back(LASteps[prevStage.End]);

		return true;
	}

	void HarmonicLLAEvaluator::Evaluate(IRasterizer rasterizer) {
		real_hr x, y;
		while (rasterizer.GetPixel(x, y)) {
			complex dc = { real(x), real(y) };
			complex Z = 0.0, z = 0.0, dz = 0.0;
			real norm_dc = chebyshev_norm(dc);

			uint_iter i = 0, j = 0;
			size_t stage = LAStages.size();
			if (LAStages.size()) j = LAStages.back().Begin;

			while (stage) {
				stage--;
				size_t begin = LAStages[stage].Begin;
				size_t end = LAStages[stage].End;

				while (i < parameters.Iterations) {
					LAStep &step = LASteps[j];
					complex new_dz = dz * (step.Z + z);

					if (chebyshev_norm(new_dz) > step.ValidRadius || norm_dc > step.ValidRadiusC) {
						j = step.NextStageLAIndex;
						break;
					}

					dz = new_dz * step.A + dc * step.B;

					i += step.Length;
					j++;

					z = dz + LASteps[j].Z;

					if (j == end || chebyshev_norm(z) < chebyshev_norm(dz)) {
						j = begin;
						dz = z;
					}
				}
			}

			Z = reference[j];
			while (i < parameters.Iterations) {
				dz = dz * (Z + z) + dc;
				i++; j++;

				Z = reference[j];
				z = Z + dz;

				if (norm(z) > 4096.0) break;

				if (j == referenceLength || norm(z) < norm(dz)) {
					Z = real(0.0);
					dz = z;
					j = 0;
				}
			}

			Output output;
			output.Value = i;

			rasterizer.WriteResults(&output);
		}
	}
}