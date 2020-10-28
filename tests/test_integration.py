from plaque_assay import main


def test_main():
    results = main.main()
    assert isinstance(results, dict)
