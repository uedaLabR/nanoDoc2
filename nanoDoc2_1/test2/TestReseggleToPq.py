from click.testing import CliRunner
import nanoDoc2_1.cmd.nanoDoc as nanoDoc

def test_reseggle():

  runner = CliRunner()
  result = runner.invoke(nanoDoc.fast5ToReSegmentedPq, ['Peter'])
  assert result.exit_code == 0
  assert result.output == 'Hello Peter!\n'