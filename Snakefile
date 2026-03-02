configfile: "config.default.yaml"
SAMPLES = list(config["samples"].keys())
import os, sys
from pprint import pformat

bamDir = config.get("bam_directory", "") or ''
if bamDir and os.path.isdir(bamDir):
    for fn in os.listdir(bamDir):
        if fn.endswith(".bam"):
            sample = fn.rsplit('.', 1)[0]
            config["samples"][sample] = os.path.join(bamDir, fn)

if not config.get("referenceGenome"):
    config["referenceGenome"] = str(workflow.basedir) + "/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
config["callerInfixes"] = [x.lower() for x in config.get("callers", [])]

print(file=sys.stderr)
print(pformat(config["samples"], sort_dicts=False), file=sys.stderr)
print(file=sys.stderr)

if not (bamDir and os.path.isdir(bamDir)):
    include: "rules/alignment.smk"

include: "rules/variants.smk"
include: "rules/insertions.smk"
include: "rules/deletions.smk"
include: "rules/reference.smk"
include: "rules/r.smk"
include: "rules/methylation.smk"
include: "rules/plotmethylation.smk"
include: "rules/met_compare.smk"
workdir: config['outputPath']
config["minimumReadSupport"] = int(config.get("minimumReadSupport", 0))

def init():
    result = [
        f"variants/{config['allPrefix']}.te.vcf.gz",
        f"variants/{config['allPrefix']}.te.lax.vcf.gz",
    ]

    # Se añaden los archivos del expand correctamente a la lista
    result.extend(
        expand("reports/meth_report.{sample}.html", sample=SAMPLES)
    )

    if config.get("mode") == "full":

        for x, y in config.get("comparisons", []):
            result.append(f"reports/{x}_vs_{y}.html")
            result.append(f"UnsharedReport/Compare_Meth_Report{x}_vs_{y}.html")

        result.append(f"reports/report.{config['allPrefix']}.html")

    return result


rule all:
    input:
        init()
