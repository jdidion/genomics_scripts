import os
from ete2a1 import PhyloTree
from ete2a1.treeview import TreeStyle

sys.path.append("%s/src/lib/python" % os.environ['LAB_HOME'])
import util.io

def make_style(style=None, mode='c', title=None, topo=False, border=False, branchlen=False, 
        branchsup=False, span=None, start=None):
    if style is None:
        style = TreeStyle()
    style.show_leaf_name = True
    style.mode = mode
    style.title = title
    style.force_topology = topo
    style.show_border = border
    style.show_branch_length = branchlen
    style.show_branch_support = branchsup
    if mode = 'c' and span is not None:
        style.arc_start = start or (-1 * span)
        style.arc_span = span

def draw_tree(newick, output=None, subformat=0, layout='phylogeny', tree_style=None, node_style=None):
    if os.path.exists(newick):
        newick = util.io.safe_read_file(newick)
    
    assert newick is not None, "Invalid Newick tree"
    
    tree = PhyloTree(newick, format=subformat)
    
    if tree_style is None:
        tree_style = make_style()
    
    if node_style is not None:
        for n in tree.traverse():
            n.set_style(node_style)
    
    if output is None:
        # draw to the screen
        tree.show(layout, tree_style)
    else:
        # write to a file
        tree.render(output, layout, tree_style=tree_style)