#include "PTWithCompression"
#include <Imagina/output_info_helper>
#include <Imagina/pixel_management>

IMPLEMENT_INTERFACE(PTWithCompression::PTWithCompressionEvaluator, Imagina::IEvaluator);

namespace PTWithCompression {
	const PixelDataInfo *PTWithCompressionEvaluator::GetOutputInfo() {
		IM_GET_OUTPUT_INFO_IMPL(Output, Value);
	}

	void PTWithCompressionEvaluator::Precompute() {
		HPComplex C = HPComplex(x, y);
		HPComplex Z = C;
		complex z = complex(Z), dzdc = 1.0;

		ReferenceCompressor compressor(reference, complex(C));

		size_t i = 1;
		while (i < parameters.Iterations) {
			compressor.Add(z);

			dzdc = 2.0 * z * dzdc + 1.0;
			Z = Z * Z + C;
			z = complex(Z);
			i++;

			if (radius * chebyshev_norm(dzdc) * 2.0 > chebyshev_norm(z)) break;
			if (norm(z) > 16.0) break;
		};

		compressor.Finalize(z);
	}

	void PTWithCompressionEvaluator::Evaluate(IRasterizingInterface rasterizingInterface) {
		HRReal x, y;
		while (rasterizingInterface.GetPixel(x, y)) {
			complex dc = { real(x), real(y) };
			complex Z = 0.0, z = 0.0, dz = 0.0;

			ReferenceDecompressor decompressor(reference);

			ITUInt i = 0;
			while (i < parameters.Iterations) {
				dz = dz * (Z + z) + dc;
				i++;

				Z = decompressor.Next();
				z = Z + dz;

				if (norm(z) > 4096.0) break;

				if (decompressor.End() || norm(z) < norm(dz)) {
					Z = decompressor.Reset();
					dz = z;
				}
			}

			Output output;
			output.Value = i;

			rasterizingInterface.WriteResults(&output);
		}
	}
}