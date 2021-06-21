import awkward

def dummy_jes_syst(events):
    """
    Dummy function illustrating a jet energy scale uncertainty that results in new jet collections with Jet.pt varied.
    Should be deleted once real examples are implemented.
    """
    jets = events.Jet 

    variations = {}
    variations["up"] = jets.pt + (10 * awkward.ones_like(jets.pt))
    variations["down"] = jets.pt - (10 * awkward.ones_like(jets.pt))
    return variations
