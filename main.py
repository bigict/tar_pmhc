"""Tools for preprocess TCR_pMHC data"""
import os
import re
from collections import defaultdict
import csv
from dataclasses import dataclass
import functools
import json
import multiprocessing as mp
import pickle
import random
import time
from urllib.parse import urlparse, parse_qs

from profold2.data.parsers import parse_fasta
from profold2.data.utils import decompose_pid, seq_index_join, seq_index_split
from profold2.utils import timing


@dataclass
class DBUri:
  path: str
  chain_idx: str = "chain.idx"
  mapping_idx: str = "mapping.idx"
  attr_idx: str = "attr.idx"
  a3m_dir: str = "a3m"


def parse_db_uri(db_uri):
  o = urlparse(db_uri)
  chain_idx, mapping_idx, attr_idx = "chain.idx", "mapping.idx", "attr.idx"
  a3m_dir = "a3m"
  if o.query:
    attrs = parse_qs(o.query)
    if "chain_idx" in attrs:
      chain_idx = attrs["chain_idx"][-1]
    if "mapping_idx" in attrs:
      mapping_idx = attrs["mapping_idx"][-1]
    if "attr_idx" in attrs:
      attr_idx = attrs["attr_idx"][-1]
    if "a3m_dir" in attrs:
      a3m_dir = attrs["a3m_dir"][-1]
  return DBUri(o.path, chain_idx, mapping_idx, attr_idx=attr_idx, a3m_dir=a3m_dir)

def _db_uri_abs_path(data_uri, data_idx):
  if os.path.isabs(data_idx):
    return data_idx
  return os.path.join(data_uri.path, data_idx)

def read_mapping_idx(data_uri, mapping_dict=None):
  if mapping_dict is None:
    mapping_dict = {}
  # mapping_idx_path = os.path.join(data_uri.path, data_uri.mapping_idx)
  mapping_idx_path = _db_uri_abs_path(data_uri, data_uri.mapping_idx)
  if os.path.exists(mapping_idx_path):
    with open(mapping_idx_path, "r") as f:
      for line in filter(lambda x: x, map(lambda x: x.strip(), f)):
        k, v = line.split()
        mapping_dict[v] = k
  return mapping_dict


def read_chain_idx(data_uri, chain_dict=None):
  if chain_dict is None:
    chain_dict = {}
  # chain_idx_path = os.path.join(data_uri.path, data_uri.chain_idx)
  chain_idx_path = _db_uri_abs_path(data_uri, data_uri.chain_idx)
  if os.path.exists(chain_idx_path):
    with open(chain_idx_path, "r") as f:
      for line in filter(lambda x: x, map(lambda x: x.strip(), f)):
        k, *v = line.split()
        chain_dict[k] = v
  return chain_dict


def read_attrs_idx(data_uri, attr_dict=None):
  if attr_dict is None:
    attr_dict = {}
  # attr_idx_path = os.path.join(data_uri.path, data_uri.attr_idx)
  attr_idx_path = _db_uri_abs_path(data_uri, data_uri.attr_idx)
  if os.path.exists(attr_idx_path):
    with open(attr_idx_path, "r") as f:
      for line in filter(lambda x: x, map(lambda x: x.strip(), f)):
        k, *v = line.split()
        try:
          attr_dict[k] = json.loads(" ".join(v))
        except: pass
  return attr_dict


def _read_fasta(data_dir, pid):
  with open(os.path.join(data_dir, "fasta", f"{pid}.fasta"), "r") as f:
    fasta_str = f.read()
  sequences, _ = parse_fasta(fasta_str)
  assert len(sequences) == 1
  return pid, sequences[0]


def read_fasta(data_uri, mapping_idx):
  fasta_list = {}
  with mp.Pool(processes=int(os.environ.get("NUM_PROCESSES", 16))) as p:
    f = functools.partial(_read_fasta, data_uri.path)
    for pid, sequence in p.imap(f,
                                set(mapping_idx.values()),
                                chunksize=int(os.environ.get("CHUNKSIZE",
                                                             256))):
      fasta_list[sequence] = pid
  return fasta_list


def create_shared_obj(**kwargs):
  context = mp.get_context()
  manager = context.Manager()
  shared_obj = manager.Namespace()
  for k, v in kwargs.items():
    setattr(shared_obj, k, v)
  return shared_obj


def align_peptide_main(args):  # pylint: disable=redefined-outer-name
  assert args.db
  seq_db, desc_db = [], []
  for db in args.db:
    with open(db, "r") as f:
      fasta_str = f.read()
    seqs, descs = parse_fasta(fasta_str)
    seq_db += seqs
    desc_db += descs
  desc_db = [desc.split()[0] for desc in desc_db]
  assert len(seq_db) == len(desc_db)
  seq_len = defaultdict(list)
  for i, seq in enumerate(seq_db):
    seq_len[len(seq)].append(i)

  for fasta_file in args.files:
    if args.verbose:
      print(f"process {fasta_file} ...")
    with open(fasta_file, "r") as f:
      fasta_string = f.read()
    sequences, descriptions = parse_fasta(fasta_string)
    assert len(sequences) == 1
    assert len(sequences) == len(descriptions)

    pid, _ = os.path.splitext(os.path.basename(fasta_file))
    output_path = os.path.join(args.output, pid, "msas")
    os.makedirs(output_path, exist_ok=True)
    with open(os.path.join(output_path, "bfd_uniclust_hits.a3m"), "w") as f:
      f.write(f">{descriptions[0]}\n")
      f.write(f"{sequences[0]}\n")

      k = len(sequences[0])
      if k in seq_len:
        for i in seq_len[k]:
          f.write(f">{desc_db[i]}/1-{k}\n")
          f.write(f"{seq_db[i]}\n")


def align_peptide_add_argument(parser):  # pylint: disable=redefined-outer-name
  parser.add_argument("files", type=str, nargs="*", help="list of fasta files")
  parser.add_argument("-o",
                      "--output",
                      type=str,
                      default=".",
                      help="output dir, default=\".\"")
  parser.add_argument("--db",
                      type=str,
                      default=None,
                      nargs="+",
                      help="peptide fasta db, default=None")
  return parser


def read_a3m(data_uri, mapping_idx, pdb_id):
  pdb_id = mapping_idx.get(pdb_id, pdb_id)
  a3m_path = os.path.join(_db_uri_abs_path(data_uri, data_uri.a3m_dir), pdb_id, "msas", f"{pdb_id}.a3m")
  with open(a3m_path, "r") as f:
    a3m_string = f.read()
  sequences, descriptions = parse_fasta(a3m_string)
  return sequences, descriptions


def align_a3m(a3m_data, mapping_dict, align_dict, **kwargs):
  a3m_seqs, a3m_descs = a3m_data

  align_data_dict, align_chain_dict = align_dict
  for seq, desc in zip(a3m_seqs[1:], a3m_descs[1:]):
    pid, chain, domains = decompose_pid(desc.split()[0], return_domain=True)
    if domains:
      domains = list(seq_index_split(domains))
    else:
      domains = []

    pdb_id = f"{pid}_{chain}" if chain else pid
    align_data_dict[pdb_id] = (domains, seq, desc, kwargs)

    pdb_id_list = set([pdb_id]) | set(mapping_dict.get(pdb_id, []))
    for pdb_id in pdb_id_list:
      pid, chain = decompose_pid(pdb_id)  # pylint: disable=unbalanced-tuple-unpacking
      # align_dict[pid].append((chain, domains, seq, desc, kwargs))
      align_chain_dict[pid].append(chain)
  # return align_dict
  return align_data_dict, align_chain_dict


def align_complex(shared_obj, target_uri, target_mapping_idx, item):
  target_pid, target_chain_list = item

  db_attr_idx = shared_obj.db_attr_idx
  db_chain_idx = shared_obj.db_chain_idx
  db_mapping_idx = shared_obj.db_mapping_idx
  db_mapping_dict = shared_obj.db_mapping_dict

  def _seq_at_i(a3m_data, i):
    seqs, descs = a3m_data
    return seqs[i], descs[i]

  # retrieve aligned chains
  a3m_list, a3m_dict = [], ({}, defaultdict(list))
  with timing(f"read_a3m_dict ({target_pid})", print_fn=print):
    for chain in target_chain_list:
      pid = f"{target_pid}_{chain}" if chain else target_pid
      a3m_data = read_a3m(target_uri, target_mapping_idx, pid)
      with timing(f"align_a3m ({target_pid}_{chain})", print_fn=print):
        a3m_dict = align_a3m(a3m_data,
                             db_mapping_dict,
                             a3m_dict,
                             target_chain=chain)
      a3m_list.append(a3m_data)

  def _is_aligned(pid, chain_list):
    if pid in db_attr_idx:
      return pid in db_chain_idx and len(db_chain_idx[pid]) == len(chain_list)
    return False

  align_data_dict, a3m_dict = a3m_dict
  # filter complex with all chains aligned
  with timing(f"filter a3m_dict ({target_pid})", print_fn=print):
    new_a3m_dict = {k: v for k, v in a3m_dict.items() if _is_aligned(k, v)}

  # realign the complex: iterate each target chain
  new_a3m_list = []

  # add target
  target_seq, domains = "", []
  n = 1
  for i, chain in enumerate(target_chain_list):
    seq, _ = _seq_at_i(a3m_list[i], 0)
    target_seq += seq
    domains.append((n, n + len(seq) - 1))
    n += len(seq) + 100
  domains = seq_index_join(domains)

  target_desc = f">{target_pid} domains:{domains}"
  new_a3m_list.append(target_desc)
  new_a3m_list.append(target_seq)

  def _repl(m):
    return "*" * len(m.group(0))

  # hit chains
  for pid, chain_list in new_a3m_dict.items():
    hit_desc = f">{pid} chains:{','.join(c for c, *_ in chain_list)}"
    if pid in db_attr_idx:
      hit_desc = f"{hit_desc} {db_attr_idx[pid]}"
    new_a3m_list.append(hit_desc)

    new_hit_seq = ""
    for i, target_chain in enumerate(target_chain_list):
      seq, _ = _seq_at_i(a3m_list[i], 0)
      hit_seq_at_i = "*" * len(seq)
      # for _, _, hit_seq, _, attrs in chain_list:
      for chain in chain_list:
        pdb_id = f"{pid}_{chain}" if chain else pid
        assert pdb_id in db_mapping_idx, (target_pid, pdb_id)
        if pdb_id in align_data_dict:
          _, hit_seq, _, attrs = align_data_dict[pdb_id]
        else:
          assert db_mapping_idx[pdb_id] in align_data_dict
          _, hit_seq, _, attrs = align_data_dict[db_mapping_idx[pdb_id]]
        if attrs["target_chain"] == target_chain:
          hit_seq_at_i = hit_seq
          hit_seq_at_i = re.sub("^[-]+", _repl, hit_seq_at_i)
          hit_seq_at_i = re.sub("[-]+$", _repl, hit_seq_at_i)
          break
      new_hit_seq += hit_seq_at_i
    new_a3m_list.append(new_hit_seq)

  return item, new_a3m_list, a3m_dict


def align_complex_main(args):  # pylint: disable=redefined-outer-name
  db_mapping_idx, db_chain_idx, db_attr_idx = {}, {}, {}
  for db_uri in args.db_uri:
    db_uri = parse_db_uri(db_uri)
    db_mapping_idx = read_mapping_idx(db_uri, db_mapping_idx)
    db_chain_idx = read_chain_idx(db_uri, db_chain_idx)
    db_attr_idx = read_attrs_idx(db_uri, db_attr_idx)
  db_mapping_dict = defaultdict(list)
  for k, v in db_mapping_idx.items():
    db_mapping_dict[v].append(k)

  target_uri = parse_db_uri(args.target_uri)
  target_mapping_idx = read_mapping_idx(target_uri)
  target_chain_idx = read_chain_idx(target_uri)

  with mp.Pool(processes=args.processes) as p:
    shared_obj = create_shared_obj(db_mapping_idx=db_mapping_idx,
                                   db_chain_idx=db_chain_idx,
                                   db_attr_idx=db_attr_idx,
                                   db_mapping_dict=db_mapping_dict)
    f = functools.partial(align_complex,
                          shared_obj,
                          target_uri,
                          target_mapping_idx)
    for (pid, chain_list), a3m_list, a3m_dict in p.imap(
        f, target_chain_idx.items(), chunksize=args.chunksize):
      output_path = os.path.join(args.output, pid, "msas")
      os.makedirs(output_path, exist_ok=True)
      with open(os.path.join(output_path, f"{pid}.a3m"), "w") as f:
        f.write("\n".join(a3m_list))
      with open(os.path.join(output_path, f"{pid}.pkl"), "wb") as f:
        pickle.dump(a3m_dict, f)
      print(f"{pid}\t{len(a3m_list)/2}\t{chain_list}")


def align_complex_add_argument(parser):  # pylint: disable=redefined-outer-name
  parser.add_argument("files", type=str, nargs="*", help="list of fasta files")
  parser.add_argument("-o",
                      "--output",
                      type=str,
                      default=".",
                      help="output dir, default=\".\"")
  parser.add_argument("--db_uri",
                      type=str,
                      default=None,
                      nargs="+",
                      help="db uri.")
  parser.add_argument("--target_uri",
                      type=str,
                      default=None,
                      help="target uri.")
  parser.add_argument("--processes",
                      type=int,
                      default=None,
                      help="num of processes.")
  parser.add_argument("--chunksize", type=int, default=2, help="chunksize.")
  return parser


def csv_to_fasta_main(args):  # pylint: disable=redefined-outer-name
  output_uri = parse_db_uri(args.output)
  os.makedirs(output_uri.path, exist_ok=True)

  print(f"load {args.target_uri} ...")
  target_uri = parse_db_uri(args.target_uri)
  mapping_idx = read_mapping_idx(target_uri)
  attr_idx = read_attrs_idx(target_uri)
  fasta_idx = read_fasta(target_uri, mapping_idx)

  def cell_check(c):
    return c != "" and c.find("nan") == -1

  def cell_write(c, pid):
    if c not in fasta_idx:
      fasta_path = os.path.join(output_uri.path, "fasta")
      os.makedirs(fasta_path, exist_ok=True)
      with open(os.path.join(fasta_path, f"{pid}.fasta"), "w") as f:
        f.write(f">{pid}\n")
        f.write(c)
      mapping_idx[pid] = pid
      fasta_idx[c] = pid
    else:
      mapping_idx[pid] = fasta_idx[c]

  for csv_file in args.csv_file:
    print(f"process {csv_file} ...")
    with open(csv_file, "r") as f:
      reader = csv.DictReader(f)

      for i, row in enumerate(reader, start=args.start_idx):
        pdb_id = f"{args.pid_prefix}{i}"

        for key, chain in (("Antigen", "P"), ("MHC_str", "M"), ("a_seq", "A"),
                           ("b_seq", "B"), ("tcrb", "B"), ("TCRA", "A"),
                           ("TCRB", "B")):
          if key in row:
            if cell_check(row[key]):
              cell_write(row[key], f"{pdb_id}_{chain}")

        if "y" in row:
          attr_idx[pdb_id] = {"label": float(row["y"])}
        elif "label" in row:
          attr_idx[pdb_id] = {"label": float(row["label"])}
        elif args.default_y is not None:
          attr_idx[pdb_id] = {"label": args.default_y}
        if "HLA" in row:
          if pdb_id in attr_idx:
            attr_idx[pdb_id]["MHC"] = row["HLA"]
          else:
            attr_idx[pdb_id] = {"MHC": row["HLA"]}

  print(f"write {output_uri.mapping_idx} ...")
  with open(os.path.join(output_uri.path, output_uri.mapping_idx), "w") as f:
    for v, k in mapping_idx.items():
      f.write(f"{k}\t{v}\n")
  print(f"write {output_uri.attr_idx} ...")
  with open(os.path.join(output_uri.path, output_uri.attr_idx), "w") as f:
    for k, v in attr_idx.items():
      v = json.dumps(v)
      f.write(f"{k}\t{v}\n")


def csv_to_fasta_add_argument(parser):  # pylint: disable=redefined-outer-name
  parser.add_argument("-o",
                      "--output",
                      type=str,
                      default=".",
                      help="output dir.")
  parser.add_argument("--target_uri", type=str, default=".", help="target dir.")
  parser.add_argument("--start_idx",
                      type=int,
                      default=0,
                      help="start index for each protein.")
  parser.add_argument("--pid_prefix",
                      type=str,
                      default="tcr_pmhc_",
                      help="pid prefix.")
  parser.add_argument("--default_y",
                      type=float,
                      default=None,
                      help="default label.")
  parser.add_argument("csv_file", type=str, nargs="+", default=None, help="csv file")
  return parser


def create_negative_main(args):  # pylint: disable=redefined-outer-name
  random.seed()

  mhc_seq_dict = {}

  peptides = set()
  tcr_mhc_rows = defaultdict(set)
  print(f"process {args.csv_file} ...")
  with open(args.csv_file, "r") as f:
    reader = csv.DictReader(f)
    for row in reader:
      tcr_mhc_rows[(row["a_seq"], row["b_seq"], row["MHC_str"])].add(
          row["Antigen"])
      peptides.add(row["Antigen"])

  print(f"write {args.output} ...")
  with open(args.output, "w") as f:
    writer = csv.DictWriter(f, ("Antigen", "a_seq", "b_seq", "MHC_str"))
    writer.writeheader()
    for i, (row, antigens) in enumerate(tcr_mhc_rows.items()):
      a_seq, b_seq, mhc_str = row
      if args.verbose:
        print(f"{i}\t{len(antigens)}/{len(peptides)}")
      negatives = list(peptides - antigens)
      if negatives:
        random.shuffle(negatives)
        for antigen in negatives[:max(1, int(len(antigens) * args.amplify))]:
          writer.writerow({
              "Antigen": antigen,
              "a_seq": a_seq,
              "b_seq": b_seq,
              "MHC_str": mhc_str
          })


def create_negative_add_argument(parser):  # pylint: disable=redefined-outer-name
  parser.add_argument("-o",
                      "--output",
                      type=str,
                      default=None,
                      help="output dir.")
  parser.add_argument("-n",
                      "--amplify",
                      type=float,
                      default=1.0,
                      help="amplify.")
  parser.add_argument("csv_file", type=str, default=None, help="csv file")
  return parser


def mhc_filter_main(args):  # pylint: disable=redefined-outer-name
  def _is_aligned(target, seq):
    i, j = 0, 0
    state = 0
    while i < len(target) and j < len(seq):
      if seq[j].islower():
        return False
      if state == 0:
        if seq[j] == "-":
          i, j = i + 1, j + 1
        elif target[i] == seq[j]:
          i, j = i + 1, j + 1
          state = 1
        else:
          return False
      elif state == 1:
        if seq[j] == "-":
          i, j = i + 1, j + 1
          state = 2
        elif target[i] == seq[j]:
          i, j = i + 1, j + 1
        else:
          return False
      elif state == 2:
        if seq[j] != "-":
          return False
        i, j = i + 1, j + 1
    # seq = seq.strip("i-")
    # if re.match(".*[-a-z].*", seq):
    #   return False
    return True

  for mhc_a3m_file in args.mhc_a3m_file:
    print(f"processing {mhc_a3m_file} ...")
    with open(mhc_a3m_file, "r") as f:
      a3m_string = f.read()
    sequences, descriptions = parse_fasta(a3m_string)
    assert len(sequences) > 0
    assert len(sequences) == len(descriptions)
    print(f"{mhc_a3m_file}\t{len(sequences)}")
    data = filter(lambda x: _is_aligned(sequences[0], x[0]),
                  zip(sequences, descriptions))
    sequences, descriptions = zip(*data)
    print(f"{mhc_a3m_file}\t{len(sequences)}")


def mhc_filter_add_argument(parser):  # pylint: disable=redefined-outer-name
  parser.add_argument("mhc_a3m_file",
                      type=str,
                      nargs="+",
                      help="mhc a3m files.")
  return parser


def mhc_preprocess_main(args):  # pylint: disable=redefined-outer-name
  mhc_seq_dict = {}

  print(f"load {args.mhc_seq_file} ...")
  with open(args.mhc_seq_file, "r") as f:
    reader = csv.DictReader(f)
    for row in reader:
      name = re.split(r"[*:]", row["name"])

      for idx in (3, 2):
        k = "".join(name[:idx])
        if k not in mhc_seq_dict:
          mhc_seq_dict[k] = row["sqe"]

  mhc_rows = []
  print(f"process {args.csv_file} ...")
  with open(args.csv_file, "r") as f:
    reader = csv.DictReader(f)

    for row in reader:
      allele = row["Allele"]
      if allele in mhc_seq_dict:
        mhc_rows.append({
            "Antigen": row["Peptide"],
            "MHC_str": mhc_seq_dict[allele]
        })
      elif args.verbose:
        print(f"{allele} not found")

  print("write {args.output} ...")
  with open(args.output, "w") as f:
    writer = csv.DictWriter(f, ("Antigen", "a_seq", "b_seq", "MHC_str"))
    writer.writeheader()
    for row in mhc_rows:
      writer.writerow(row)


def mhc_preprocess_add_argument(parser):  # pylint: disable=redefined-outer-name
  parser.add_argument("-o",
                      "--output",
                      type=str,
                      default=".",
                      help="output dir.")
  parser.add_argument("--mhc_seq_file",
                      type=str,
                      default=None,
                      help="mhc sequence file.")
  parser.add_argument("csv_file", type=str, default=None, help="csv file")
  return parser


def split_data_main(args):  # pylint: disable=redefined-outer-name
  random.seed()

  def mapping_filter(k):
    _, chain = decompose_pid(k)  # pylint: disable=unbalanced-tuple-unpacking
    return chain == args.cluster_chain

  print(f"load {args.target_uri} ...")
  target_uri = parse_db_uri(args.target_uri)
  mapping_idx = read_mapping_idx(target_uri)
  mapping_idx = {k: v for k, v in mapping_idx.items() if mapping_filter(k)}
  fasta_idx = read_fasta(target_uri, mapping_idx)
  fasta_idx = {v: k for k, v in fasta_idx.items()}  # pid -> fasta
  chain_idx = read_chain_idx(target_uri)

  cluster_idx = {}
  with open(args.cluster_csv_file, "r") as f:
    for row in csv.DictReader(f):
      cluster_idx[row["Antigen"]] = row["Cluster"]

  # split by cluster id
  cluster_list = list(set(cluster_idx.values()))
  random.shuffle(cluster_list)
  print(f"Total clusters: {len(cluster_list)}")
  test_cluster_list = set(
      cluster_list[:int(len(cluster_list) * args.test_ratio)])
  print(f"Clusters for test: {len(test_cluster_list)}")
  with open(os.path.join(target_uri.path, "cluster.idx"), "w") as f:
    for seq, c in cluster_idx.items():
      if c in test_cluster_list:
        f.write(f"{seq},{c},test\n")
      else:
        f.write(f"{seq},{c},train\n")

  for pid, chain_list in chain_idx.items():
    label = "train"
    if args.cluster_chain in chain_list:
      k = f"{pid}_{args.cluster_chain}"
      k = mapping_idx.get(k, k)
      c = fasta_idx[k]
      if c in cluster_idx and cluster_idx[c] in test_cluster_list:
        label = f"test\t{cluster_idx[c]}"
    print(f"split_data\t{pid} {' '.join(chain_list)}\t{label}")


def split_data_add_argument(parser):  # pylint: disable=redefined-outer-name
  parser.add_argument("--target_uri", type=str, default=".", help="target uri.")
  parser.add_argument("--test_ratio",
                      type=float,
                      default=0.2,
                      help="test set ratio.")
  parser.add_argument("--cluster_chain",
                      type=str,
                      default="P",
                      help="cluster by this chain only.")
  parser.add_argument("cluster_csv_file",
                      type=str,
                      default=None,
                      help="cluster csv file")
  return parser


def tcr_pmhc_to_pmhc_main(args):  # pylint: disable=redefined-outer-name
  output_uri = parse_db_uri(args.output)
  os.makedirs(output_uri.path, exist_ok=True)

  print(f"load {args.target_uri} ...")
  target_uri = parse_db_uri(args.target_uri)

  chain_idx = read_chain_idx(target_uri)
  mapping_idx = read_mapping_idx(target_uri)
  attr_idx = read_attrs_idx(target_uri)

  new_mapping_dict = defaultdict(list)
  for pid, chain_list in chain_idx.items():
    if pid in attr_idx and attr_idx[pid]["label"] > 0:
      if "P" in chain_list and "M" in chain_list:
        new_mapping_dict[(mapping_idx[f"{pid}_P"], mapping_idx[f"{pid}_M"])].append(pid)

  for (pid_p, pid_m), new_pid_list in new_mapping_dict.items():
    print(f"new_pid_list_count: {pid_p} {pid_m} {len(new_pid_list)}")

    m = len(new_pid_list)
    if 0 < args.pid_topk < m:
      m = args.pid_topk

    for i, new_pid in enumerate(new_pid_list):
      weight = m / len(new_pid_list)
      attr_idx[new_pid].update(weight=weight)

      if len(chain_idx[new_pid]) >= 3:
        new_pid = f"{args.pid_prefix}{new_pid}"

        mapping_idx[f"{new_pid}_P"] = pid_p
        mapping_idx[f"{new_pid}_M"] = pid_m

        attr_idx[new_pid] = {"label":1.0, "weight":weight}

  print(f"write {output_uri.mapping_idx} ...")
  with open(os.path.join(output_uri.path, output_uri.mapping_idx), "w") as f:
    for v, k in mapping_idx.items():
      f.write(f"{k}\t{v}\n")
  print(f"write {output_uri.attr_idx} ...")
  with open(os.path.join(output_uri.path, output_uri.attr_idx), "w") as f:
    for k, v in attr_idx.items():
      v = json.dumps(v)
      f.write(f"{k}\t{v}\n")


def tcr_pmhc_to_pmhc_add_argument(parser):  # pylint: disable=redefined-outer-name
  parser.add_argument("-o",
                      "--output",
                      type=str,
                      default=".",
                      help="output dir.")
  parser.add_argument("--target_uri", type=str, default=".", help="target dir.")
  parser.add_argument("--start_idx",
                      type=int,
                      default=0,
                      help="start index for each protein.")
  parser.add_argument("--pid_prefix",
                      type=str,
                      default="pmhc_",
                      help="pid prefix.")
  parser.add_argument("--pid_topk",
                      type=int,
                      default=1,
                      help="pid prefix.")
  return parser


if __name__ == "__main__":
  import argparse

  commands = {
      "align_peptide": (align_peptide_main, align_peptide_add_argument),
      "align_complex": (align_complex_main, align_complex_add_argument),
      "csv_to_fasta": (csv_to_fasta_main, csv_to_fasta_add_argument),
      "create_negatives": (create_negative_main, create_negative_add_argument),
      "mhc_preprocess": (mhc_preprocess_main, mhc_preprocess_add_argument),
      "mhc_a3m_filter": (mhc_filter_main, mhc_filter_add_argument),
      "split_data": (split_data_main, split_data_add_argument),
      "tcr_pmhc_to_pmhc": (tcr_pmhc_to_pmhc_main, tcr_pmhc_to_pmhc_add_argument),
  }

  formatter_class = argparse.ArgumentDefaultsHelpFormatter
  parser = argparse.ArgumentParser(formatter_class=formatter_class)

  sub_parsers = parser.add_subparsers(dest="command", required=True)
  for cmd, (_, add_argument) in commands.items():
    cmd_parser = sub_parsers.add_parser(cmd, formatter_class=formatter_class)
    cmd_parser = add_argument(cmd_parser)
    cmd_parser.add_argument("-v",
                            "--verbose",
                            action="store_true",
                            help="verbose")

  args = parser.parse_args()
  main, _ = commands[args.command]
  main(args)
