import sys

class Node(object):

    def __init__(self, blen=0.0, name="", left=None, right=None):
        self.blen = blen
        self.name = name
        self.left = left
        self.right = right

    def __repr__(self):
        if len(self.name) == 0:
            return "internode:%.1g" % self.blen
        else:
            return "%s:%.1g" % (self.name, self.blen)

class Counter(object):
    def __init__(self):
        self.count = 0

    def inc(self):
        self.count += 1
        return self.count

def read_tree(nwk):
    name = ""
    blen = ""
    left = None
    right = None
    read_blen = False
    stack = [Node()]
    idx = Counter()
    while idx.count < len(nwk):
        c = nwk[idx.count]
        if c == ';':
            break
        if c in ('\n \t'):   # ignore whitespace
            idx.inc()
            continue
        if c == '(': # read left node
            idx.inc()
            stack.append(Node())
        elif c == ')': # done reading branch length
            read_blen = False
            idx.inc()
            stack[-1].blen = float(blen)
            blen = ""
            tmp = stack.pop()
            stack[-1].right = tmp
            #break
        elif c == ',':
            read_blen = False
            stack[-1].blen = float(blen)
            blen = ""
            tmp = stack.pop()
            stack[-1].left = tmp
            stack.append(Node())
            idx.inc()
        elif c == ':':  # read branch length
            if len(name) > 0:
                stack[-1].name = name
                name = ""
            read_blen = True
            idx.inc()
        else:           # reading branch lenght or name
            if read_blen:
                blen += c
            else:
                name += c
            idx.inc()
    return stack[-1]

def read_nwk(nwk_path):
    with open(nwk_path, 'r') as f:
        nwk_str = f.read()
    return read_tree(nwk_str)



def collapse(node):
    ret = None
    if not (node.left is None or node.right is None):
        ret = node
    elif not (node.left is None and node.right is None):
        tmp = node.right if node.left is None else node.left
        tmp.blen += node.blen
        ret = tmp
    return ret

def descend_left(node, parents, side):
    tmp = node
    while tmp is not None:
        parents.append(tmp)
        tmp = tmp.left
        if tmp:
            side.append('l')    # tmp is the left child of its parent

def print_node(node):
    if len(node.name):
        print "node -  %s:%.1g" % (node.name, node.blen)
    else:
        print "node -  internode:%.1g" % (node.blen)

def trim_tree(root, leaf_valid):
    parents = [root]
    side = ['l']
    descend_left(root.left, parents, side)
    while side:
        node = parents.pop()
        if node is None:
            continue
        if len(node.name) > 0:
            tmp = node if leaf_valid(node) else None
            if side[-1] == 'l':
                parents[-1].left = tmp
            else:
                parents[-1].right = tmp
        else:
            child = side.pop()
            if child == 'l':
                parents.append(node)
                side.append('r')
                descend_left(node.right, parents, side)
            else:
                tmp = collapse(node)
                if len(side) == 0:
                    break
                if side[-1] == 'l':
                    parents[-1].left = tmp
                else:
                    parents[-1].right = tmp

    root = collapse(root)
    root.blen = 0.0
    return root

def _nwk_descend_left(node, stack, nwk):
    prev = None
    while node is not None:
        stack.append(node)
        prev = node
        node = node.left
        if node is not None:
            nwk += "("
    nwk += "%s:%.5f" % (prev.name, prev.blen)
    stack.pop()
    return nwk

def to_nwk(node):
    ret = ""
    parents = list()
    ret = _nwk_descend_left(node, parents, ret)

    last_split = list()
    #last_split.append(node)
    while parents:
        tmp = parents.pop()
        if len(last_split) > 0 and last_split[-1] == tmp:
            ret += "):%.5f" % tmp.blen
            last_split.pop()
            continue
        if len(tmp.name) == 0: # internode
            ret += ","
            last_split.append(tmp)
            parents.append(tmp)
            ret = _nwk_descend_left(tmp.right, parents, ret)
        #else:
            #if len(parents) > 0 and parents[-1] == tmp.right:
            #    ret += "):%.2g" % tmp.blen
    ret += ";"
    return ret

def to_str(node):
    if node is None:
        return 'None'
    if len(node.name) > 0:
        return "%s:%.2g" % (node.name, node.blen)
    else:
        lstr = to_str(node.left)
        rstr = to_str(node.right)
        return "(%s,%s):%.2g" % (lstr, rstr, node.blen)

if __name__ == '__main__':

    def get_filter(s):
        def func(oid):
            return oid.name[0] in s
        return func

    def test(nwk, fstr):
        root = read_tree(nwk)
        print "keeping %s from %s" % (fstr, to_str(root))
        root = trim_tree(root, get_filter(fstr))
        print "result: %s %s" % (to_str(root), to_nwk(root))
        return root

    with open(sys.argv[1], 'r') as f:
        nwk_str = f.read()

    node = test(nwk_str, 'DCB')
    test(nwk_str, 'CB')
    test(nwk_str, 'AD')
    test(nwk_str, 'CD')

    print to_nwk(node)
