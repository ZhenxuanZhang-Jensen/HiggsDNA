def lhe_scale0_sf(events, year, central_only, input_collection, working_point = "none"):
    variations = {}
    variations["central"] = events.LHEScaleWeight_Four

    if not central_only:
        variations["up"] 		= events.LHEScaleWeight_Five
        variations["down"] 	= events.LHEScaleWeight_Three

    return variations


def lhe_scale1_sf(events, year, central_only, input_collection, working_point = "none"):
    variations = {}
    variations["central"] = events.LHEScaleWeight_Four

    if not central_only:
        variations["up"] 		= events.LHEScaleWeight_Seven
        variations["down"] 	= events.LHEScaleWeight_One

    return variations

def lhe_scale2_sf(events, year, central_only, input_collection, working_point = "none"):
    variations = {}
    variations["central"] = events.LHEScaleWeight_Four

    if not central_only:
        variations["up"] 		= events.LHEScaleWeight_Eight
        variations["down"] 	= events.LHEScaleWeight_Zero

    return variations

def lhe_pdf_sf(events, year, central_only, input_collection, working_point = "none"):
    variations = {}
    variations["central"] = events.LHEPdfWeight_Unit

    if not central_only:
        variations["up"] 		= events.LHEPdfWeight_Up
        variations["down"] 	= events.LHEPdfWeight_Down

    return variations
