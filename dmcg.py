#!/usr/bin/env python3

import os, glob
os.environ['PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION'] = 'python'
import copy
import os,sys,gzip
import argparse
from rdkit.Chem import AllChem as Chem
import torch
import torch.nn.functional as F
from torch_geometric.data import Data, InMemoryDataset

import numpy as np
import random
from confgen.model.gnn import GNN
from confgen.molecule.graph import rdk2graph
from confgen.utils.utils import  set_rdmol_positions
import io, subprocess


parser = argparse.ArgumentParser()
parser.add_argument("--smi", type=str, required=True,help="Input SMILES")
parser.add_argument("--out",type=str, required=True,help="Output SDF file name")
parser.add_argument("--maxconfs",type=int,default=1,help="Number of conformers to generate")
parser.add_argument("--device", type=int, default=0)
parser.add_argument("--global-reducer", type=str, default="sum",help=argparse.SUPPRESS)
parser.add_argument("--node-reducer", type=str, default="sum",help=argparse.SUPPRESS)
parser.add_argument("--graph-pooling", type=str, default="sum",help=argparse.SUPPRESS)
parser.add_argument("--dropedge-rate", type=float, default=0.1,help=argparse.SUPPRESS)
parser.add_argument("--dropnode-rate", type=float, default=0.1,help=argparse.SUPPRESS)
parser.add_argument("--num-layers", type=int, default=6,help=argparse.SUPPRESS)
parser.add_argument("--decoder-layers", type=int, default=None,help=argparse.SUPPRESS)
parser.add_argument("--latent-size", type=int, default=256,help=argparse.SUPPRESS)
parser.add_argument("--mlp-hidden-size", type=int, default=1024,help=argparse.SUPPRESS)
parser.add_argument("--mlp_layers", type=int, default=2,help=argparse.SUPPRESS)
parser.add_argument("--use-layer-norm", action="store_true", default=False,help=argparse.SUPPRESS)

parser.add_argument("--batch-size", type=int, default=128,help=argparse.SUPPRESS)
parser.add_argument("--epochs", type=int, default=100,help=argparse.SUPPRESS)
parser.add_argument("--num-workers", type=int, default=0,help=argparse.SUPPRESS)
# parser.add_argument("--log-dir", type=str, default="", help="tensorboard log directory")
parser.add_argument("--checkpoint-dir", type=str, default="",help=argparse.SUPPRESS)

parser.add_argument("--log-interval", type=int, default=100,help=argparse.SUPPRESS)
parser.add_argument("--dropout", type=float, default=0.1,help=argparse.SUPPRESS)
parser.add_argument("--encoder-dropout", type=float, default=0.0,help=argparse.SUPPRESS)
parser.add_argument("--lr", type=float, default=1e-4,help=argparse.SUPPRESS)
parser.add_argument("--layernorm-before", action="store_true", default=False,help=argparse.SUPPRESS)
parser.add_argument("--use-bn", action="store_true", default=True,help=argparse.SUPPRESS)
parser.add_argument("--weight-decay", type=float, default=1e-2,help=argparse.SUPPRESS)
parser.add_argument("--use-adamw", action="store_true", default=True,help=argparse.SUPPRESS)
parser.add_argument("--beta2", type=float, default=0.999,help=argparse.SUPPRESS)
parser.add_argument("--period", type=float, default=10,help=argparse.SUPPRESS)

parser.add_argument("--base-path", type=str, default="")
parser.add_argument(
    "--dataset-name", type=str, default="drugs", choices=["qm9", "drugs", "iso17"],help=argparse.SUPPRESS
)
parser.add_argument("--train-size", type=float, default=0.8,help=argparse.SUPPRESS)
parser.add_argument("--seed", type=int, default=2021)
parser.add_argument("--lr-warmup", action="store_true", default=True,help=argparse.SUPPRESS)
parser.add_argument("--enable-tb", action="store_true", default=False,help=argparse.SUPPRESS)
parser.add_argument("--aux-loss", type=float, default=0.0,help=argparse.SUPPRESS)
parser.add_argument("--train-subset", action="store_true", default=True,help=argparse.SUPPRESS)
parser.add_argument("--eval-from", type=str, default=None)
parser.add_argument(
    "--data-split", type=str, choices=["cgcf", "default", "confgf"], default="default",help=argparse.SUPPRESS)
parser.add_argument("--reuse-prior", action="store_true", default=True,help=argparse.SUPPRESS)
parser.add_argument("--cycle", type=int, default=1,help=argparse.SUPPRESS)

parser.add_argument("--vae-beta", type=float, default=1.0,help=argparse.SUPPRESS)
parser.add_argument("--eval-one", action="store_true", default=False,help=argparse.SUPPRESS)
parser.add_argument("--workers", type=int, default=20,help=argparse.SUPPRESS)
parser.add_argument("--extend-edge", action="store_true", default=False,help=argparse.SUPPRESS)
parser.add_argument("--use-ff", action="store_true", default=False,help=argparse.SUPPRESS)
parser.add_argument("--threshold", type=float, default=0.5,help=argparse.SUPPRESS)
parser.add_argument("--pred-pos-residual", action="store_true", default=True,help=argparse.SUPPRESS)
parser.add_argument("--node-attn", action="store_true", default=True,help=argparse.SUPPRESS)
parser.add_argument("--global-attn", action="store_true", default=False,help=argparse.SUPPRESS)
parser.add_argument("--shared-decoder", action="store_true", default=False,help=argparse.SUPPRESS)
parser.add_argument("--shared-output", action="store_true", default=True,help=argparse.SUPPRESS)
parser.add_argument("--sample-beta", type=float, default=1.2,help=argparse.SUPPRESS)
parser.add_argument("--remove-hs", action="store_true", default=True,help=argparse.SUPPRESS)
parser.add_argument("--prop-pred", action="store_true", default=False,help=argparse.SUPPRESS)

parser.add_argument("--score", action="store_true", default=False,help=argparse.SUPPRESS)
parser.add_argument("--sigma-begin", type=float, default=10.0,help=argparse.SUPPRESS)
parser.add_argument("--sigma-end", type=float, default=0.01,help=argparse.SUPPRESS)
parser.add_argument("--noise-level", type=int, default=10,help=argparse.SUPPRESS)
parser.add_argument("--noise-steps", type=int, default=100,help=argparse.SUPPRESS)
parser.add_argument("--noise-lr", type=float, default=2.4e-6,help=argparse.SUPPRESS)
parser.add_argument("--decoder-std", type=float, default=1.0,help=argparse.SUPPRESS)
parser.add_argument("--score-prior", action="store_true", default=False,help=argparse.SUPPRESS)

args = parser.parse_args()

np.random.seed(args.seed)
torch.manual_seed(args.seed)
torch.cuda.manual_seed(args.seed)
random.seed(args.seed)

device = (
    torch.device("cuda:" + str(args.device))
    if torch.cuda.is_available()
    else torch.device("cpu")
)
shared_params = {
    "mlp_hidden_size": args.mlp_hidden_size,
    "mlp_layers": args.mlp_layers,
    "latent_size": args.latent_size,
    "use_layer_norm": args.use_layer_norm,
    "num_message_passing_steps": args.num_layers,
    "global_reducer": args.global_reducer,
    "node_reducer": args.node_reducer,
    "dropedge_rate": args.dropedge_rate,
    "dropnode_rate": args.dropnode_rate,
    "dropout": args.dropout,
    "layernorm_before": args.layernorm_before,
    "encoder_dropout": args.encoder_dropout,
    "use_bn": args.use_bn,
    "vae_beta": args.vae_beta,
    "decoder_layers": args.decoder_layers,
    "reuse_prior": args.reuse_prior,
    "cycle": args.cycle,
    "pred_pos_residual": args.pred_pos_residual,
    "node_attn": args.node_attn,
    "global_attn": args.global_attn,
    "shared_decoder": args.shared_decoder,
    "sample_beta": args.sample_beta,
    "shared_output": args.shared_output,
}
model = GNN(**shared_params).to(device)

if not args.eval_from:
    args.eval_from = os.path.join(args.base_path,"Large_Drugs/checkpoint_94.pt")
checkpoint = torch.load(args.eval_from, map_location=device)["model_state_dict"]
cur_state_dict = model.state_dict()
del_keys = []
for k in checkpoint.keys():
    if k not in cur_state_dict:
        del_keys.append(k)
for k in del_keys:
    del checkpoint[k]
model.load_state_dict(checkpoint)
model.eval()

if args.out.endswith('.gz'):
    out = gzip.open(args.out,'wt')
else:
    out = open(args.out,'wt')
writer = Chem.SDWriter(out)
confs = args.maxconfs

with open(args.smi,'rt') as f:
    for smi in f:
        print(smi)
        try:
            mol = Chem.MolFromSmiles(smi)
            mol = Chem.AddHs(mol)
            Chem.EmbedMultipleConfs(mol)
            if mol.GetNumConformers() == 0:  #failed, try w/o sc
                Chem.RemoveStereochemistry(mol)
                Chem.EmbedMultipleConfs(mol) # there is exactly one platinum molecule where this triggers
            mol = Chem.RemoveAllHs(mol)
            data = Data()
            graph = rdk2graph(mol)
            data.edge_index = torch.from_numpy(graph["edge_index"]).to(torch.int64)
            data.edge_attr = torch.from_numpy(graph["edge_attr"]).to(torch.int64)
            data.x = torch.from_numpy(graph["node_feat"]).to(torch.int64)
            data.n_nodes = graph["n_nodes"]
            data.n_edges = graph["n_edges"]
            data.pos = torch.from_numpy(mol.GetConformer(0).GetPositions()).to(torch.float)

            data.rd_mol = copy.deepcopy(mol)
            data.nei_src_index = torch.from_numpy(graph["nei_src_index"]).to(torch.int64)
            data.nei_tgt_index = torch.from_numpy(graph["nei_tgt_index"]).to(torch.int64)
            data.nei_tgt_mask = torch.from_numpy(graph["nei_tgt_mask"]).to(torch.bool)
            data.num_graphs = 1
            batch, slices = InMemoryDataset.collate([data])

            for _ in range(confs):
                with torch.no_grad():
                    atoms, _ = model(batch.cuda(),sample=True)
                mol = set_rdmol_positions(data.rd_mol,atoms[-1])
                writer.write(mol)
        except Exception as e:
            print("Error",e,"with",smi,flush=True)

writer.close()
out.close() 

if args.out.endswith('.sdf.gz'):
    base = args.smi.replace('.smi','')
    refsdf = base.replace('_nosc','')+'.sdf'
    out = args.out.replace('.sdf.gz','.rmsds.txt')
    text = subprocess.check_output(f'obrms -f -m {refsdf} {args.out} > {out}',shell=True)

