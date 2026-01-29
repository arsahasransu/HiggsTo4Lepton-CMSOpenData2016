# Combine json file per dataset to a single json file

import json


def combine_json_files(newfilename: str, oldfilelist: list[dict]):

    newevents = []
    duplicate_event_count = 0

    for oldfiledict in oldfilelist:

        oldfiletype = oldfiledict['type']
        oldfile = oldfiledict['file']

        with open(oldfile, 'r') as f:
            events = json.load(f)

        for event in events:
            event['type'] = oldfiletype
            if len(newevents) == 0:
                newevents.append(event)
            else:
                found_event = False
                for ex_events in newevents:
                    if event['event'] == ex_events['event']:
                        print(event)
                        found_event = True
                        duplicate_event_count += 1
                if not found_event:
                        # print(f"Append new event to exisitng size {len(newevents)}")
                        newevents.append(event)

    with open(f"../skimmingcrabconfigs/{newfilename}", 'w') as f:
        print(f"Count of duplicate events: {duplicate_event_count}")
        print(f"Dumping merged json of {len(newevents)} events "
              f"to ../skimmingcrabconfigs/{newfilename}")
        json.dump(newevents, f, indent = 4)


if __name__ == "__main__":

    oldfilelist = [{'type': '4mu', 'file': '4mu_doublemu_2016G.json'},
                   {'type': '2mu2e', 'file': '2mu2e_doublemu_2016G.json'}]
    combine_json_files("doublemu_2016G.json", oldfilelist)

    oldfilelist = [{'type': '4mu', 'file': '4mu_doublemu_2016H.json'},
                   {'type': '2mu2e', 'file': '2mu2e_doublemu_2016H.json'}]
    combine_json_files("doublemu_2016H.json", oldfilelist)

    oldfilelist = [{'type': '4mu', 'file': '4mu_singlemu_2016G.json'},
                   {'type': '2mu2e', 'file': '2mu2e_singlemu_2016G.json'}]
    combine_json_files("singlemu_2016G.json", oldfilelist)

    oldfilelist = [{'type': '4mu', 'file': '4mu_singlemu_2016H.json'},
                   {'type': '2mu2e', 'file': '2mu2e_singlemu_2016H.json'}]
    combine_json_files("singlemu_2016H.json", oldfilelist)

    oldfilelist = [{'type': '4e', 'file': '4e_doubleel_2016G.json'},
                   {'type': '2mu2e', 'file': '2mu2e_doubleel_2016G.json'}]
    combine_json_files("doubleel_2016G.json", oldfilelist)

    oldfilelist = [{'type': '4e', 'file': '4e_doubleel_2016H.json'},
                   {'type': '2mu2e', 'file': '2mu2e_doubleel_2016H.json'}]
    combine_json_files("doubleel_2016H.json", oldfilelist)

    oldfilelist = [{'type': '4e', 'file': '4e_singleel_2016G.json'},
                   {'type': '2mu2e', 'file': '2mu2e_singleel_2016G.json'}]
    combine_json_files("singleel_2016G.json", oldfilelist)

    oldfilelist = [{'type': '4e', 'file': '4e_singleel_2016H.json'},
                   {'type': '2mu2e', 'file': '2mu2e_singleel_2016H.json'}]
    combine_json_files("singleel_2016H.json", oldfilelist)

    oldfilelist = [{'type': '2mu2e', 'file': '2mu2e_mueg_2016G.json'}]
    combine_json_files("mueg_2016G.json", oldfilelist)

    oldfilelist = [{'type': '2mu2e', 'file': '2mu2e_mueg_2016H.json'}]
    combine_json_files("mueg_2016H.json", oldfilelist)
