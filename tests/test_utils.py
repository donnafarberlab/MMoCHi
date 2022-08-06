import pytest
from mmochi import utils

def test_list_tools():
    assert utils.list_tools([1,2,3],'+',[4,5,6]) == [1,2,3,4,5,6]
