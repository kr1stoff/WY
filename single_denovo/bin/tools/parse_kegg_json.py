import json
import re
import sys
import os

if len(sys.argv) != 2:
    print(f"python3 parse_json.py <kegg.json>")
    exit(-1)

injson = os.path.abspath(sys.argv[1])
out = os.path.dirname(injson)
with open(injson) as f:
    ko_map_data = json.load(f)

with open(f"{out}/KEGG_pathway_ko.txt", "w") as oh:
    line = "#level1_pathway_id\tlevel1_pathway_name\tlevel2_pathway_id\tlevel2_pathway_name"
    line += "\tlevel3_pathway_id\tlevel3_pathway_name\tko\tko_name\tko_des\tec\n"
    oh.write(line)
    for level1 in ko_map_data["children"]:
        m = re.match(r"(\S+)\s+([\S\w\s]+)", level1["name"])
        level1_pathway_id = m.groups()[0].strip()
        level1_pathway_name = m.groups()[1].strip()
        for level2 in level1["children"]:
            m = re.match(r"(\S+)\s+([\S\w\s]+)", level2["name"])
            level2_pathway_id = m.groups()[0].strip()
            level2_pathway_name = m.groups()[1].strip()
            for level3 in level2["children"]:
                #m = re.match(r"(\S+)\s+([^\[]*)", level3["name"])  #不取'['后的信息
                m = re.match(r"(\S+)\s+(.*)", level3["name"])
                level3_pathway_id = m.groups()[0].strip()
                level3_pathway_name = m.groups()[1].strip()
                if "children" in level3:
                    for ko in level3["children"]:
                        m = re.match(r"(\S+)\s+(\S+);\s+([^\[]+)\s*(\[EC:\S+(?:\s+[^\[\]]+)*\])*", ko["name"])
                        if m is not None:
                            ko_id = m.groups()[0].strip()
                            ko_name = m.groups()[1].strip()
                            ko_des = m.groups()[2].strip()
                            ec = m.groups()[3]
                            if ec==None:
                                ec = "-"
                        line = level1_pathway_id + "\t" + level1_pathway_name + "\t" + level2_pathway_id + "\t" + level2_pathway_name
                        line += "\t" + level3_pathway_id + "\t" + level3_pathway_name + "\t" + ko_id + "\t" + ko_name + "\t" + ko_des + "\t" + ec + "\n"
                        oh.write(line)