import pandas as pd
from urllib.parse import unquote
from collections import defaultdict

def parse_attributes(attributes):
    key_value_pairs = tuple(pair_str.split('=', 1) for pair_str in attributes.split(';'))
    return {key_value[0]: unquote(key_value[1]) for key_value in key_value_pairs}

def get_descendants(parent_to_children, parent):
    descendants = []
    for child in parent_to_children[parent]:
        descendants.append(child)
        descendants += get_descendants(parent_to_children, child)
    descendants_without_duplicates = list(dict.fromkeys(descendants))
    return descendants_without_duplicates

def get_ancestors(child_to_parents, child):
    ancestors = []
    for parent in child_to_parents[child]:
        ancestors.append(parent)
        ancestors += get_ancestors(child_to_parents, parent)
    ancestors_without_duplicates = list(dict.fromkeys(ancestors))
    return ancestors_without_duplicates

def get_parent_to_descendants(parent_to_children):
    return {parent: get_descendants(parent_to_children, parent) for parent in parent_to_children}

def get_child_to_ancestors(child_to_parents):
    return {child: get_ancestors(child_to_parents, child) for child in child_to_parents}
    
def get_parent_to_children(child_to_parents):
    parent_to_children = {ID: [] for ID in child_to_parents}
    for ID, parents in child_to_parents.items():
        for parent in parents:
            parent_to_children[parent].append(ID)
    return parent_to_children

def read_gff(path, parse_attrs=True, infer_relations=True, identifier='ID'):
    gff = pd.read_csv(path, sep='\t', comment='#', header=None)
    gff.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    gff['length'] = gff.end - gff.start + 1
    if parse_attrs:
        attributes_table = pd.DataFrame(parse_attributes(attributes) for attributes in gff.attributes)
        gff = pd.concat([gff, attributes_table], axis=1)
        gff.drop(columns=['attributes'], inplace=True)
        if infer_relations and 'Parent' in gff.columns:
            child_to_parents = {
                ID: parents.split(',') if not pd.isna(parents) else []
                for ID, parents in zip(gff[identifier], gff.Parent)
            }
            child_to_ancestors = get_child_to_ancestors(child_to_parents)
            parent_to_children = get_parent_to_children(child_to_parents)
            parent_to_descendants = get_parent_to_descendants(parent_to_children)
            gff['Parents'] = [child_to_parents[ID] for ID in gff.ID]
            gff['Ancestors'] = [child_to_ancestors[ID] for ID in gff.ID]
            gff['Children'] = [parent_to_children[ID] for ID in gff.ID]
            gff['Descendants'] = [parent_to_descendants[ID] for ID in gff[identifier]]
    return gff
