import os
import logging

import torch
from torch.utils.data.distributed import DistributedSampler

from profold2.command import worker
from profold2.data import dataset
from profold2.model.head import HeaderBuilder


def _load_models(rank, args):
  def _location_split(model_location):
    k = model_location.find("=")
    if k != -1:
      return model_location.split("=", 1)
    model_name = os.path.basename(model_location)
    model_name, _ = os.path.splitext(model_name)
    return model_name, model_location

  for i, model_location in enumerate(args.refs):
    model_name, model_location = _location_split(model_location)
    logging.info("Load model [%d/%d] %s from %s",
        i, len(refs), model_name, model_location)
    pth = torch.load(model_location)
    yield model_name, pth

def _create_dataloader(rank, args):
  kwargs = {"pin_memory": True, "shuffle": False}
  if xpu.is_available() and WorkerXPU.world_size(args.nnodes) > 1:
    kwargs["num_replicas"] = WorkerXPU.world_size(args.nnodes)
    kwargs["rank"] = xpu.rank
  return dataset.load(
      data_dir=args.data_dir,
      data_idx=args.data_idx,
      chain_idx=f"{args.vars}/chain.idx",
      mapping_idx=f"{args.vars}/mapping.idx",
      var_dir=f"{args.vars}/var",
      pseudo_linker_prob=1.0,
      pseudo_linker_shuffle=False,
      num_workers=args.num_workers, **kwargs)

def predict(rank, args):
  dim_single, dim_pairwise = 384, 128
  headers = [("fitness", {"mask": "-", "pooling": "mean"}, {})]

  headers = HeaderBuilder.build((dim_single, dim_pairwise), headers)

  # model_runners = dict(_load_models(rank, args))

  test_loader = _create_dataloader(rank, args)

  def predict_fitness(idx, batch):
    pass

  for idx, batch in enumerate(iter(test_loader)):
    try:
      predict_fitness(idx, batch)
    except RuntimeError as e:
      logging.error("%d %s", idx, str(e))
    


def add_arguments(parser):  # pylint: disable=redefined-outer-name
  cwd = os.path.dirname(__file__)

  parser.add_argument("--refs", type=str, nargs="+",
                      help="list of reference pth files")
  parser.add_argument("--db", type=str,
                      default=os.path.join(cwd, "data/tcr_pmhc_db"),
                      help="reference data dir")
  parser.add_argument("--vars", type=str, default=None,
                      help="aligned ratio threshold.")


if __name__ == "__main__":
  import argparse

  formatter_class = argparse.ArgumentDefaultsHelpFormatter
  parser = argparse.ArgumentParser(formatter_class=formatter_class)

  # init distributed env
  parser.add_argument("--nnodes", type=int, default=None,
      help="number of nodes.")
  parser.add_argument("--node_rank", type=int, default=0,
      help="rank of the node.")
  parser.add_argument("--local_rank", type=int, default=None,
      help="local rank of xpu.")
  parser.add_argument("--init_method", type=str,
      default="file:///tmp/profold2-tcr_pmhc.dist",
      help="method to initialize the process group")

  # output dir
  parser.add_argument("-o", "--prefix", type=str, default=".",
      help="prefix of out directory.")
  add_arguments(parser)
  # verbose
  parser.add_argument("-v", "--verbose", action="store_true", help="verbose")

  args = parser.parse_args()
  worker.main(args, predict)
