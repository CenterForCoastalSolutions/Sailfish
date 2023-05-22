# collin to make a draft of mod sources. then connect with jose again (just outline the structure,functions etc, doesnt need to be perfect)


#for some reason in fortran (ROMS) they set a maximum number of sources manually to 200. No clue what it's about
#it seems like it's just a check in case a user put a value of like 2000 sources by accident but there is no explaination in the code.
ANA_PSOURCE_2D = True
#for now, only implementing analytical sources, can do netcdf later



def apply_sources(SOURCES):
    #seperate 2D and 3D as a 3D source could be in the intermediate layers of the water
    if ANA_PSOURCE_2D:
        SOURCES.Isrc

