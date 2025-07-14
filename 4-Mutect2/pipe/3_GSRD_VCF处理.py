#!/usr/bin/env python3
import argparse
import pysam

def main():
    """
    主函数：处理多等位点VCF文件，仅保留支持度最高的ALT等位基因，强制样本基因型为1/1，并只输出GT字段。

    功能说明：
    1. 从输入VCF文件中，对于每个位点，如果存在多个ALT等位基因，选择AD支持度最高的ALT。
    2. 将样本的基因型（GT）强制设为1/1。
    3. 删除除GT以外的所有FORMAT字段，仅保留GT。
    4. 输出处理后的VCF文件。

    参数:
        -i, --input: 输入VCF文件路径（必需）
        -o, --output: 输出VCF文件路径（必需）

    注意事项：
    - 假设VCF文件中仅包含一个样本。
    - 依赖pysam库进行VCF文件的读取与写入。
    """
    parser = argparse.ArgumentParser(
        description="从多等位 VCF 中保留支持度最高的 ALT，强制 GT=1/1，并只输出 GT 字段"
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="输入 VCF 文件路径"
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="输出 VCF 文件路径"
    )
    args = parser.parse_args()

    in_vcf_path = args.input
    out_vcf_path = args.output

    # 打开输入／输出
    vcf_in  = pysam.VariantFile(in_vcf_path,  "r")
    vcf_out = pysam.VariantFile(out_vcf_path, "w", header=vcf_in.header)

    # 假设只有一个样本
    sample = list(vcf_in.header.samples)[0]

    for rec in vcf_in:
        # 1. 如果多等位，选支持度最高的那个 ALT
        ad = rec.samples[sample]['AD']
        ref_count = ad[0]
        alt_counts = ad[1:]
        # best_index 是 alt_counts 中最大值的索引
        best_index = max(range(len(alt_counts)), key=lambda i: alt_counts[i])
        rec.alts = (rec.alts[best_index],)
        best_alt_count = alt_counts[best_index]

        # 2. 强制 GT=1/1
        rec.samples[sample]['GT'] = (1, 1)

        # 3. 删除除 GT 之外的所有 FORMAT 字段
        for key in list(rec.samples[sample].keys()):
            if key != 'GT':
                del rec.samples[sample][key]

        # 写出
        vcf_out.write(rec)

    vcf_in.close()
    vcf_out.close()
    print("Done →", out_vcf_path)

if __name__ == "__main__":
    main()

# python3 filter_biallelic.py -i /path/to/input.vcf -o /path/to/output.vcf
