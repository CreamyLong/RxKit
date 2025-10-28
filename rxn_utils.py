from mol_utils import canonicalize_smiles, inchify, get_longest_smiles


def process_reaction(rxn, with_agent=True):
    """
    Process and canonicalize reaction SMILES
    """
    reactants, reagents, products = rxn.split(">")

    if with_agent:
        try:
            precursors = [inchify(canonicalize_smiles(r)) for r in reactants.split(".")]
            if len(reagents) > 0:
                precursors += [
                    inchify(canonicalize_smiles(r)) for r in reagents.split(".")
                ]
            products = [inchify(canonicalize_smiles(p)) for p in products.split(".")]
            precursors_set = set(precursors)
            products_set = set(products)
            common_set = precursors_set & products_set
            products_set = products_set - common_set
            precursors = list(precursors_set)
            products = list(products_set)

        except Exception as E:
            # print("process_reaction",E)
            return ""

    elif not with_agent:
        try:
            precursors = [inchify(canonicalize_smiles(r)) for r in reactants.split(".")]
            products = [inchify(canonicalize_smiles(p)) for p in products.split(".")]
        except Exception as E:
            return ""

    joined_precursors = ".".join(sorted(precursors, key=lambda s: len(s), reverse=True))
    joined_products = get_longest_smiles(".".join(sorted(products, key=lambda s: len(s), reverse=True)))
    return f"{joined_precursors}>>{joined_products}"
