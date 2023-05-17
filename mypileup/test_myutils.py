import argparse
import os
import sys
import pytest

from . import myutils

# TODO - add unit tests here

def test_GetQual():
	assert(myutils.GetQual(0)=="!")
