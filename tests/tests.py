import unittest
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../cobreeder_for_ghost_wolves')))
import ghost_cobreeder_v2_CLI as cli


class TestArgparse(unittest.TestCase):
    def test_main_normal(self, mock_print):
        cli.main()
        mock_print.assert_called_once_with("Hello, Charlie!")

    def test_main_shout(self, mock_print):
        cli.main()
        mock_print.assert_called_once_with("HELLO, DANA!")

if __name__ == "__main__":
    unittest.main()