import jgidb
import sys
from trim_tree import trim_tree, read_nwk

nwk_path = sys.argv[1]
root = read_nwk(nwk_path)

oids = list()
def get_leaf_name(node):
    oids.append(node.name)
    return True

trim_tree(root, get_leaf_name)

def get_query(oids):
    sql = "SELECT t.taxon_oid, t.jgi_project_id FROM taxon t WHERE t.taxon_oid in (%s)"
    sub = ",".join(str(x) for x in oids)
    return sql % sub

step = 100
for i in range(0, len(oids), step):
    query = get_query(oids[i:i+step])
    for oid, spid in jgidb.queryDb("IMG", query):
        if spid in (0, None):
            continue
        print "%d\t%s" % (oid, str(spid))

