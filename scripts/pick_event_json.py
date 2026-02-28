# Select events using a json mask
# Assumes json format produced by DumpJson module
# Produces a merged art file <dataset>.pick_events.art

import os
import json
from argparse import ArgumentParser

#-------------------------------------------------------------------------------------
# Merge the individual event files
#-------------------------------------------------------------------------------------
def merge_events(dataset):
    base = '.'.join(dataset.split('.')[0:-1])
    out_file = base + '.pick_events' + '.' + dataset.split('.')[-1]
    file_list = f'{dataset}.pick_event.files'
    file_pattern = dataset + '_*_*_*.art'
    command = f'ls {file_pattern} >| {file_list}'
    print(command)
    os.system(command)
    command = f'art -c /dev/null -o {out_file} -S {file_list}'
    print(command)
    os.system(command)

#-------------------------------------------------------------------------------------
# Retrieve the event as an art file
#-------------------------------------------------------------------------------------
def pick_event(run, subrun, event, dataset):
    command = f'pickEvent -e {dataset} {run}/{subrun}/{event}'
    print(command)
    os.system(command
    return

#-------------------------------------------------------------------------------------
# Main function
#-------------------------------------------------------------------------------------
if __name__ == "__main__":

    # Get the inputs
    parser = ArgumentParser()
    parser.add_argument(
        "-j",
        "--json",
        dest="json",
        default="skim_trig_clusters.json",
        help="Input json mask file",
    )
    parser.add_argument(
        "-d",
        "--dataset",
        dest="dataset",
        required=True,
        help="Input dataset",
    )
    parser.add_argument(
        "-m",
        "--max-events",
        dest="max_events",
        default=-1,
        type=int,
        help="Maximum events to process",
    )

    args = parser.parse_args()

    json_file  = args.json
    dataset    = args.dataset
    max_events = args.max_events

    print(f"Using input dataset {dataset} with json mask file {json_file}")

    mask_file = open(json_file, 'r')
    mask = json.load(mask_file)

    # Process the events from the mask
    nevents = 0
    for run_entry in mask:
        run = run_entry['run']
        for subrun_entry in run_entry['subRuns']:
            subrun = subrun_entry['subRun']
            for event in subrun_entry['events']:
                pick_event(run, subrun, event, dataset)
                nevents += 1
                if max_events > 0 and nevents >= max_events: break
            if max_events > 0 and nevents >= max_events: break
        if max_events > 0 and nevents >= max_events: break

    # Merge the events
    merge_events(dataset)
