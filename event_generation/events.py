base_path = "/work/jeppelt/girlsday2022/new_event_prod"


event_dict = {
    "mumu": ["13", "-13"],
    "pipi": ["211", "-211"],
    "eepipi": ["11", "-11", "211", "-211"],
    "ee": ["11", "-11"],
    "eemumu": ["11", "-11", "13", "-13"],
    "gammagamma": ["22", "22"],
    "eegamma": ["11", "-11", "22"],
    "mumugamma": ["13", "-13", "22"],
    "pipigamma": ["211", "-211", "22"],
    "DM0ee": ["d0", "11", "-11"],
    "DM0gammaee": ["d0", "22", "-11", "11"],
    "DM4DM4eemumu": ["d4", "d4", "-11", "11", "-13", "13"],
    "DM1ee": ["d1", "-11", "11"],
    "DM2ee": ["d2", "-11", "11"],
    "DM2pipi": ["d2", "-211", "211"],
    "DM3mumu": ["d3", "-13", "13"]
}
events_dict = {}

for i in range(1):
    for key, value in event_dict.items():
        events_dict.update({f"{key}{i}": value})