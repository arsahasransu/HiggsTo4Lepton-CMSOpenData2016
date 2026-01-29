import json
import matplotlib.pyplot as plt
import os
import random
import yaml
import zipfile


def postmix():

    print("Hello World")


def event_exists(event_list, event):
    
    if len(event_list) == 0:
        return False
    
    exists = any([event['event']==e['event'] for e in event_list])
    if exists:
        # print("Skipping duplicate event: \n", event)
        return True
    
    return False


def find_event_in_igfiles(event, dfolder):
    igfiles = os.listdir(dfolder)

    for igfile in igfiles:

        with zipfile.ZipFile(dfolder+igfile, "r") as z:

            for eventname in z.namelist():
                if str(event['event']) in eventname:
                    return dfolder+igfile    


def make_unique_events(sigsets):

    unique_events = []

    for sigset in sigsets:
        json_file = sigset['json']
        dfolder = sigset['igfiles']

        with open(json_file) as f:
            events = json.load(f)

        for event in events:
            if not event_exists(unique_events, event):

                # Find IGFile with event
                igfilename = find_event_in_igfiles(event, dfolder)

                # Append filename name to event
                event['file'] = igfilename
                unique_events.append(event)

    return unique_events


def make_combined_mass_plot(sevents):

    masses = []
    n4mu = 0
    n2mu2e = 0
    n4e = 0

    for e in sevents:
        masses.append(e["fourlep_mass"])
        if e['type'] == '4mu':
            n4mu += 1
        elif e['type'] == '2mu2e':
            n2mu2e += 1
        else:
            n4e += 1
    fig = plt.figure()
    
    counts, bins, _ = plt.hist(masses, 108, range=(70, 502), histtype='step')
    bin_centers = (bins[:-1] + bins[1:]) / 2

    print(f'\nRatio of signal composition: ({n4mu}:{n2mu2e}:{n4e}) (4mu:2mu2e:4e)')
    plt.plot(bin_centers, counts, 'o', label="CMS Outreach 2016, $16.6\\ fb^{-1}$")
    plt.xlabel('m(l+l-l+l-) / GeV')
    plt.ylabel('Events / 4 GeV')
    plt.xscale('log')
    plt.legend()
    plt.savefig("combined_mass_hist.png")


class BackgroundSet:
    def __init__(self, fpath):
        bkgigfilenames = os.listdir(fpath)
        self.bkgeventnames = []
        self.bkgevents = []
        self.current_event_ctr = 0

        for bkgigfilename in bkgigfilenames:
            with zipfile.ZipFile(fpath+bkgigfilename, 'r') as bkgigfile:
                bkgignames = bkgigfile.namelist()
                for bkgigname in bkgignames:
                    if bkgigname == 'Header':
                        continue

                    self.bkgeventnames.append(bkgigname)
                    with bkgigfile.open(bkgigname) as bkgigevent:
                        bkgevent = json.load(bkgigevent)
                    self.bkgevents.append(bkgevent)

    def get_bkg_len(self):
        return len(self.bkgevents)
    
    def get_next_background_event(self):
        if self.current_event_ctr < len(self.bkgeventnames):
            current_event_name = self.bkgeventnames[self.current_event_ctr]
            current_event = self.bkgevents[self.current_event_ctr]
            self.current_event_ctr += 1
        else:
            raise RuntimeError('Exceeding available background events')
        
        return (current_event_name, current_event)


def make_shuffled_ig_sets(shuffled_events, nsets, nevtpset, bkgset):
    
    for nset in range(nsets):
        shuffled_event_set = shuffled_events[nset*nevtpset:(nset+1)*nevtpset]
        
        with zipfile.ZipFile(f'./mixedigfiles/fourlepton_{nset}.ig',
                             "w", compression=zipfile.ZIP_DEFLATED) as zoutset:
            for shuffled_event in shuffled_event_set:
                inigfile = shuffled_event['file']

                if inigfile == 'background':
                    igbkgname, igbkgevent = bkgset.get_next_background_event()
                    igeventdump = json.dumps(igbkgevent)
                    zoutset.writestr(igbkgname, igeventdump)
                else:
                    run = shuffled_event['run']
                    event = shuffled_event['event']
                    name = f'Events/Run_{run}/Event_{event}'

                    with zipfile.ZipFile(inigfile, "r") as inigf:
                        with inigf.open(name) as inigeventf:
                            inigevent = json.load(inigeventf)

                    igeventdump = json.dumps(inigevent)
                    zoutset.writestr(name, igeventdump)
        
        with open(f'./mixedigfiles/event_info_{nset}.json', 'w') as json_dump_f:
            json.dump(shuffled_event_set, json_dump_f, indent=4)


if __name__ == "__main__":

    with open('dataconfig.yml') as f:
        dataconfig = yaml.safe_load(f)

    nsets = dataconfig['General']['sets']
    nevtpset = dataconfig['General']['eventsperset']
    print(nsets, nevtpset)

    sigsets = dataconfig['SignalSets']
    print(len(sigsets))

    bkgsetfpath = dataconfig['BackgroundSet']['igfiles']
    bkgset = BackgroundSet(bkgsetfpath)

    sevents = make_unique_events(sigsets)
    print("Length of unique signal events: ", len(sevents))

    make_combined_mass_plot(sevents)
    
    # Gather the necessary background events
    nsevents = len(sevents)
    nbevents = nsets*nevtpset - nsevents
    nbevents = nbevents if nbevents > 0 else 0

    print("\nRequired number of background events: ", nbevents)
    nallbevents = bkgset.get_bkg_len()
    print("Existing size of background set: ", nallbevents)
    if nbevents > nallbevents:
        raise RuntimeError(f'Required background event count {nbevents} larger than size of existing background set {nallbevents}')
    else:
        print(f'Condition met - background events {nbevents}/{nallbevents} (Required/Available)')

    bevents = [{'file': 'background'}] * nbevents

    allevents = sevents + bevents
    shuffled_events = random.sample(allevents, len(allevents))
    
    print("\nShuffling datasets and generating outreach data!\n")
    make_shuffled_ig_sets(shuffled_events, nsets, nevtpset, bkgset)
