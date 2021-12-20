import  argparse

args = argparse.ArgumentParser()
args.add_argument('--dataset', default='Caudoviridae')
args.add_argument('--model', default='gcn')
args.add_argument('--learning_rate', type=float, default=0.01)
args.add_argument('--epochs', type=int, default=200)
args.add_argument('--hidden', type=int, default=64)
args.add_argument('--dropout', type=float, default=0)
args.add_argument('--weight_decay', type=float, default=5e-4)
args.add_argument('--max_degree', type=int, default=3)
args.add_argument('--gpus', type=int, default = 0)
args.add_argument('--t', type=float, default=0.0)


args = args.parse_args()
print(args)