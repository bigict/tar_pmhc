import os
import sys
import multiprocessing as mp
import re

from profold2.data.parsers import parse_fasta, parse_a3m

def _read_a3m_name_list(filename):
  def lines(f):
    for line in filter(lambda x: len(x)>0, map(lambda x: x.strip(), f)):
      if line.startswith('>'):
        yield line[1:]

    name_list = []

    with open(filename, 'r') as f:
      n = ''
      for i, line in enumerate(lines(f)):
        m = p.search(line)
        if m:
          w = float(m.group(1))
        else:
          w = 1.0
        name_list.append((i, w, line))

    return filename, name_list

def extract_pid_main(args):  # pylint: disable=redefined-outer-name
  p = re.compile("weight[\'\"]\s*:\s*([0-9.]+)")

  if args.files:
    filename_list = args.files
  else:
    filename_list = list(filter(lambda x: len(x)>0,
                                map(lambda x: x.strip(), sys.stdin)))

  with mp.Pool() as p:
    for filename, name_list in p.imap(_read_a3m_name_list,
                                      filename_list,
                                      chunksize=4):
      for i, w, line in name_list:
        print(f'{filename}\t{i}\t{w}\t{line}')


def extract_pid_add_argument(parser):  # pylint: disable=redefined-outer-name
  parser.add_argument('files', type=str, nargs='*',
                      help='list of a3m files')
  return parser


def _aligned_ratio(msa, query=None, n=None, trim_gap=True):
  i, j = 0, len(msa)
  if trim_gap:
    while i < len(msa) and msa[i] == "-":
      i += 1
    while j > 0 and msa[j - 1] == "-":
      j -= 1
  if not n:
    n = len(msa[i:j])
  if query:
    c = 0
    while i < j:
      if query[i] != msa[i]:
        c += 1
      i += 1
  else:
    c = msa[i:j].count("-")
  r = 1. - c / n
  return r

def filter_main(args):  # pylint: disable=redefined-outer-name
  for fasta_file in args.files:
    print(f"processing {fasta_file} ...")
    pid, _ = os.path.splitext(os.path.basename(fasta_file))
    a3m_file = os.path.join(args.output, pid, 'msas', f'{pid}.a3m')
    with open(a3m_file, "r") as f:
      a3m_string = f.read()

    sequences, descriptions = parse_fasta(a3m_string)
    aligned_seqs, _ = parse_a3m(a3m_string)

    with open(a3m_file, "w") as f:
      n = 0
      for i, seq in enumerate(aligned_seqs):
        msa_aligned_ratio = _aligned_ratio(seq, query=sequences[0], trim_gap=args.trim_gap)
        if msa_aligned_ratio >= args.aligned_ratio_threshold:
          f.write(f">{descriptions[i]} aligned_ratio:{msa_aligned_ratio}\n")
          f.write(f"{sequences[i]}\n")
          n += 1
      print(f"filtering {fasta_file} {n}/{len(sequences)}")

def filter_add_argument(parser):  # pylint: disable=redefined-outer-name
  parser.add_argument('files', type=str, nargs='*',
                      help='list of fasta files')
  parser.add_argument('-o', '--output', type=str, default='.',
                      help='output dir, default=\'.\'')
  parser.add_argument('-t', '--aligned_ratio_threshold', type=float, default=0.5,
                      help='aligned ratio threshold.')
  parser.add_argument('--trim_gap', action='store_true',
                      help='trim gap.')
  return parser

if __name__ == '__main__':
  import argparse

  commands = {
      "filter": (filter_main, filter_add_argument),
      "extract_pid": (extract_pid_main, extract_pid_add_argument),
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
