import os

GLOBAL_COLLECTION = "/mnt/N/NGS/TB_pipeline/TB_pipeline_database/DB"
LOCAL_COLLECTION="/media/nacho/Data/temp/toptest/tempdb"

globaldirs = next(os.walk(GLOBAL_COLLECTION))[1]
localdirs=next(os.walk(LOCAL_COLLECTION))[1]

if len(globaldirs) == 0:
    globaldirstxt = ""
else:
    globaldirstxt = " ".join([GLOBAL_COLLECTION + "/" + g for g in globaldirs])

if len(localdirs) == 0:
    localdirstxt = ""
else:
    localdirstxt=" ".join([LOCAL_COLLECTION + "/" + l for l in localdirs])
    globaldirstxt=globaldirstxt+localdirstxt



f = open("demofile3.txt", "w")
f.write(globaldirstxt)
f.close()
