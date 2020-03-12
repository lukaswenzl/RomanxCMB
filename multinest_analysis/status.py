import argparse
import watcher

#example: python modules/WFIRSTxCMB/multinest_analysis/status.py --n "mn_quick_DESY1"

parser = argparse.ArgumentParser(description='Getting the status from multinest')
parser.add_argument("--n", default="No Name Given please specify root", type=str, help="Root name for the files of the current multinest run")
parser.add_argument("--p", default=False, type=bool, help="Option to make some plots")
parser.add_argument("--n_params", default=16, type=int, help="Number of parameters for plots")



args = parser.parse_args()

root = args.n

mn_printer = watcher.ProgressPrinter(-1, outputfiles_basename = root)
mn_printer.run_once()

if(args.p):
    mn_plotter = watcher.ProgressPlotter(args.n_params, outputfiles_basename = root)
    mn_plotter.run_once()
    print("plots saved")



