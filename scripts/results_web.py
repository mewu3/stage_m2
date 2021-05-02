#!/usr/bin/env python3.8
import sys
import re
import os
from collections import defaultdict
from Bio import SeqIO

# python3 scripts/evaluation-basic.py /datater/wu/data/MSSPE-basic/enterovirus/kmer15/allOligo.set /datater/wu/data/MSSPE-basic/enterovirus/enterovirus.msa /datater/wu/data/MSSPE-basic/enterovirus/kmer15/evaluation/report.html /datater/wu/data/MSSPE-basic/enterovirus/kmer15/evaluation/style.css enterovirus


inputOligo = sys.argv[1]
inputMSA = sys.argv[2]
htmlFile = sys.argv[3]
cssFile = sys.argv[4]
virusName = sys.argv[5]

records = list(SeqIO.parse(inputMSA, "fasta"))
seqLength = len(records[0].seq)
# virusName = "ZIKA virus"

dict_kmerInfo = defaultdict(list)

with open(inputOligo, "r") as file:
    for line in file.readlines()[1:]:
        line = line.rstrip("\n").split() # index oligo kmerCount CG% Entropy Tm homodimer-dG hairpin-dG start end chunk-start chunk-end
        dict_kmerInfo[line[0]]=line[1:]

# page HTML ####################################################################
htmlFileOpen = open(htmlFile, "w")

formattedContent = """{{
                            accesibility: {{description: "{nature}"}},
                            name: "{sequence}",
                            description: "{description}",
                            x: {start},
                            x2: {end},
                            y: {order},
                        }},"""

dataContent = ""

dataContent += formattedContent.format(
    nature = 'genome',
    sequence = records[0].id,
    description = records[0].seq,
    start = 0,
    end = seqLength,
    order = 0
)

for id in dict_kmerInfo:
    info = dict_kmerInfo[id]
    formatedList = f"kmerCount:{info[1]}, CG%:{info[2]}, Entropy:{info[3]}, Tm:{info[4]}"
    content = formattedContent.format(
        nature = 'oligonucleotide',
        sequence = info[0],
        description=formatedList,
        start = int(info[-1]) - 50 + 1,
        end = info[-1],
        order = 1
    )
    dataContent += content

htmlFileOpen.write(f"""<!doctype html>
<html>
    <head>
        <meta charset="utf-8">
        <title>
            test
        </title>
        <link rel="stylesheet" href="style.css">
        <script src="https://code.highcharts.com/highcharts.js"></script>
        <script src="https://code.highcharts.com/modules/xrange.js"></script>
        <script src="https://code.highcharts.com/modules/exporting.js"></script>
        <script src="https://code.highcharts.com/modules/accessibility.js"></script>
        <script>https://code.highcharts.com/css/highcharts.css</script>
    </head>
    <body>
        <figure class="highcharts-figure">
            <div id="container"></div>
        </figure>
        <script>
        document.addEventListener('DOMContentLoaded', function(){{
            const chart = Highcharts.chart('container', {{
                chart: {{
                    type: 'xrange',
                    styledMode: true
                }},
                title: {{
                    text: "MSSPE report"
                }},
                yAxis: {{
                    title: {{
                        text: ''
                    }},
                    categories: ['Genome', 'oligonucleotide']
                }},
                tooltip: {{
                    useHTML: true,
                    pointFormat : (
                        '{{point.name}}<br/>'+'{{point.description}}'
                    ),
                }},
                series: [{{
                    borderColor: 'gray',
                    pointWidth: 20,
                    data: [{
                        dataContent
                    }],
                    dataLabels: {{
                        enabled:true
                    }}
                }}]
            }});
        }})
        </script>
    </body>
</html>""")

htmlFileOpen.close()

# fichier CSS ##################################################################
cssFileOpen = open(cssFile, "w")

cssFileOpen.write('''@import 'https://code.highcharts.com/css/highcharts.css';
#container {
    width: 100%;
    height: 400px;
    margin: 1em auto;
}
.highcharts-xrange-series .highcharts-point {
	stroke-width: 1px;
	stroke: gray;
}
.highcharts-partfill-overlay {
	fill: #010101;
	fill-opacity: 0.3;
}
.highcharts-data-label text {
	fill: white;
	text-shadow: 1px 1px black, -1px 1px black, -1px -1px black, 1px -1px black;
}
.highcharts-tooltip {
    word-wrap: break-word;
    font-size: 10px ;
}
#genomeSeq {
    word-wrap: break-word;
    font-size: 10px;
    display: block;
    height: 300px;
    overflow: scroll;
}''')

cssFileOpen.close()

################################################################################
#!/usr/bin/env python3.8
import sys
import re
import os
from collections import defaultdict
from Bio import SeqIO

# python3 scripts/evaluation-variant.py /datater/wu/data/MSSPE-variant/enterovirus/kmer15/allOligo.set /datater/wu/data/MSSPE-variant/enterovirus/enterovirus0.99.msa /datater/wu/data/MSSPE-variant/enterovirus/kmer15/evaluation/report.html /datater/wu/data/MSSPE-variant/enterovirus/kmer15/evaluation/style.css enterovirus

# python3 scripts/evaluation-variant.py /datater/wu/data/MSSPE-variant/coronavirus/kmer15/allOligo.set /datater/wu/data/MSSPE-variant/coronavirus/coronavirus0.99.msa /datater/wu/data/MSSPE-variant/coronavirus/kmer15/evaluation/report.html /datater/wu/data/MSSPE-variant/coronavirus/kmer15/evaluation/style.css coronavirus


inputOligo = sys.argv[1]
inputMSA = sys.argv[2]
htmlFile = sys.argv[3]
cssFile = sys.argv[4]
virusName = sys.argv[5]

records = list(SeqIO.parse(inputMSA, "fasta"))
seqLength = len(records[0].seq)
# virusName = "ZIKA virus"

dict_kmerInfo = defaultdict(list)

with open(inputOligo, "r") as file:
    for line in file.readlines()[1:]:
        line = line.rstrip("\n").split() # index oligo kmerCount CG% Entropy Tm homodimer-dG hairpin-dG start end chunk-start chunk-end
        dict_kmerInfo[line[0]]=line[1:]

# page HTML ####################################################################
htmlFileOpen = open(htmlFile, "w")

formattedContent = """{{
                            accesibility: {{description: "{nature}"}},
                            name: "{sequence}",
                            description: "{description}",
                            x: {start},
                            x2: {end},
                            y: {order},
                        }},"""

dataContent = ""

dataContent += formattedContent.format(
    nature = 'genome',
    sequence = records[0].id,
    description = records[0].seq,
    start = 0,
    end = seqLength,
    order = 0
)

for id in dict_kmerInfo:
    info = dict_kmerInfo[id]
    formatedList = f"kmerCount:{info[1]}, CG%:{info[2]}, Entropy:{info[3]}, Tm:{info[4]}"
    content = formattedContent.format(
        nature = 'oligonucleotide',
        sequence = info[0],
        description=formatedList,
        start = info[-4],
        end = info[-3],
        order = 1
    )
    dataContent += content
    content2 = formattedContent.format(
        nature = 'oligonucleotide',
        sequence = info[0],
        description=formatedList,
        start = info[-2],
        end = info[-1],
        order = 2
    )
    dataContent += content2

htmlFileOpen.write(f"""<!doctype html>
<html>
    <head>
        <meta charset="utf-8">
        <title>
            test
        </title>
        <link rel="stylesheet" href="style.css">
        <script src="https://code.highcharts.com/highcharts.js"></script>
        <script src="https://code.highcharts.com/modules/xrange.js"></script>
        <script src="https://code.highcharts.com/modules/exporting.js"></script>
        <script src="https://code.highcharts.com/modules/accessibility.js"></script>
        <script>https://code.highcharts.com/css/highcharts.css</script>
    </head>
    <body>
        <figure class="highcharts-figure">
            <div id="container"></div>
        </figure>
        <script>
        document.addEventListener('DOMContentLoaded', function(){{
            const chart = Highcharts.chart('container', {{
                chart: {{
                    type: 'xrange',
                    styledMode: true
                }},
                title: {{
                    text: "MSSPE report"
                }},
                yAxis: {{
                    title: {{
                        text: ''
                    }},
                    categories: ['Genome', 'oligonucleotide-readPosition', 'oligonucleotide']
                }},
                tooltip: {{
                    useHTML: true,
                    pointFormat : (
                        '{{point.name}}<br/>'+'{{point.description}}'
                    ),
                }},
                series: [{{
                    borderColor: 'gray',
                    pointWidth: 20,
                    data: [{
                        dataContent
                    }],
                    dataLabels: {{
                        enabled:true
                    }}
                }}]
            }});
        }})
        </script>
    </body>
</html>""")

htmlFileOpen.close()

# fichier CSS ##################################################################
cssFileOpen = open(cssFile, "w")

cssFileOpen.write('''@import 'https://code.highcharts.com/css/highcharts.css';
#container {
    width: 100%;
    height: 400px;
    margin: 1em auto;
}
.highcharts-xrange-series .highcharts-point {
	stroke-width: 1px;
	stroke: gray;
}
.highcharts-partfill-overlay {
	fill: #010101;
	fill-opacity: 0.3;
}
.highcharts-data-label text {
	fill: white;
	text-shadow: 1px 1px black, -1px 1px black, -1px -1px black, 1px -1px black;
}
.highcharts-tooltip {
    word-wrap: break-word;
    font-size: 10px ;
}
#genomeSeq {
    word-wrap: break-word;
    font-size: 10px;
    display: block;
    height: 300px;
    overflow: scroll;
}''')

cssFileOpen.close()
