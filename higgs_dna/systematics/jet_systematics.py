import awkward

def dummy_jes_syst(events, is_data):
    """
    Dummy function illustrating a jet energy scale uncertainty that results in new jet collections with Jet.pt varied.
    Should be deleted once real examples are implemented.
    """
    jets = events.Jet 

    variations = {}
    variations["central"] = jets.pt + (2 * awkward.ones_like(jets.pt))
    if not is_data:
        variations["up"] = jets.pt + (12 * awkward.ones_like(jets.pt))
        variations["down"] = jets.pt - (8 * awkward.ones_like(jets.pt))
    return variations
