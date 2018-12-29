import nxpd
import pydot
import networkx as nx

def draw_circuit(circ, circ_name):

    G = nx.MultiDiGraph()
    filename = circ_name + '.png'

    for gate_id, gate in enumerate(circ.op_tab):

        g_name = []
        for qg in gate.gate_list:
            g_name.append(qg.name)    

        G.add_node(gate_id)

        G.node[gate_id]['name'] = ''.join(g_name) + str(gate_id)
    
    for gate_id, gate in enumerate(circ.op_tab):
        for wire in gate.wires:
            for g_suc in gate.suc[wire]:

                G.add_edge(gate_id, circ.op_tab.index(g_suc), name = wire)
    
    G.graph["dpi"] = 70 
    for node in G.nodes(data = True):
        n = G.nodes[node[0]]
        n['label'] = node[1]['name'] 
        n['color'] = 'blue'
        n['style'] = 'filled'
        n['fillcolor'] = 'lightblue'
    for e in G.edges(data=True):
        e[2]['label'] = e[2]['name']

    show = True 

    return nxpd.draw(G, filename=filename, show=show)
