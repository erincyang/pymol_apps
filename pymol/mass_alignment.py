def mass_cealign(target):
    list = cmd.get_names("public_objects")
    filter(lambda x:cmd.get_type(x)!="object:molecule",list)
    for name in list:
        if name!=target:
           cmd.cealign('polymer and name ca and (%s)'%target,'polymer and name ca and (%s)'%name)

def mass_tmalign(target):
    list = cmd.get_names("public_objects")
    filter(lambda x:cmd.get_type(x)!="object:molecule",list)
    for name in list:
        if 'polar_conts' in name:
            continue
        if name!=target:
           cmd.tmalign('polymer and name ca and (%s)'%name, 'polymer and name ca and (%s)'%target)

def mass_pair_fit(ligand_selection, req_lig_name=None, align_to_first=True):
    target_atom_names = {'names':[], 'object': []}
    cmd.iterate(ligand_selection, 'names.append(name)', space=target_atom_names)
    cmd.iterate(ligand_selection, 'object.append(model)', space=target_atom_names)
    sel_str = 'name ' + ' name '.join(target_atom_names['names'])

    # Make the target selection
    name = target_atom_names['object'][0]
    target_atoms = name + ' and (' + sel_str + ')'
    #target_atoms += f' and (byres first resn C2E and {name}) '
    cmd.select('tmp_target', target_atoms)
    print(target_atoms)
    
    list = cmd.get_names("public_objects")
    filter(lambda x:cmd.get_type(x)!="object:molecule",list)
    for name in list:
        if 'polar_conts' in name:
            continue
        
        # Make the movable selection
        movable_atoms = name + ' and (' + sel_str + ') ' # and chain A
        #movable_atoms += f' and (byres first resn C2E and {name})'
        print(movable_atoms)
        cmd.select('tmp', movable_atoms)
        cmd.pair_fit('tmp', 'tmp_target')


cmd.extend("mass_cealign", mass_cealign)
cmd.extend("mass_tmalign", mass_tmalign)
cmd.extend("mass_pair_fit", mass_pair_fit)

