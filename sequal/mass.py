from sequal import resources


def calculate_mass(seq, mass_dict=None, N_terminus=0, O_terminus=0, with_water=True):
    """
    Calculates the mass of a Sequence object.

    Args:
        seq (sequal.sequence.Sequence): A Sequence object representing the
            sequence of amino acids.
        mass_dict (dict, optional): A dictionary containing the masses of
            potential modifications and amino acids. Defaults to None.
        N_terminus (int or float, optional): The mass at the N-terminus of the
            sequence. Defaults to 0.
        O_terminus (int or float, optional): The mass at the C-terminus of the
            sequence. Defaults to 0.
        with_water (bool, optional): Whether or not to add the mass of water.
            Defaults to True.

    Returns:
        float: The calculated mass of the sequence.

    Raises:
        ValueError: If a block or modification mass is not available in the
            mass attribute and no additional mass_dict was supplied.
    """
    mass = 0
    if with_water:
        mass += resources.H * 2 + resources.O
    for i in seq:
        if not i.mass:
            if mass_dict:
                if i.value in mass_dict:
                    mass += mass_dict[i.value]
                else:
                    raise ValueError("Block {} not found in mass_dict".format(i.value))
            else:
                raise ValueError(
                    "Block {} mass is not available in mass attribute and no additional mass_dict was supplied".format(
                        i.value
                    )
                )
        else:
            mass += i.mass
        if i.mods:
            for m in i.mods:
                if m.mass != 0 and not m.mass:
                    if mass_dict:
                        if m.value in mass_dict:
                            mass += mass_dict[m.value]
                        else:
                            raise ValueError(
                                "Block {} not found in mass_dict".format(m.value)
                            )
                    else:
                        raise ValueError(
                            "Block {} mass is not available in mass attribute and no additional mass_dict was supplied".format(
                                m.value
                            )
                        )
                else:
                    mass += m.mass
    return mass + N_terminus + O_terminus
