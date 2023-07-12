#!/usr/bin/env python
import sys
import os
import subprocess
from myparsing import *
from myutils import *

argvals = None

def main(argvals=argvals):
    P = Parser(argvals)
    version = set_version_number()
    print(f"\ngwaswa v{version} \n")

    args = sys.argv[1:]
    args.insert(0, 'python3')
    args.insert(1, os.path.join(os.path.dirname(__file__), 'mystarter.py'))



    process = subprocess.run(args)
    return process.returncode


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt as e:
        print("\n interrupted")
    finally:
        pass