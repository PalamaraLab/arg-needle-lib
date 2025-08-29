import arg_needle_lib
import pytest

def test_adna_volume():
    arg = arg_needle_lib.ARG(0, 100, 6)
    arg.add_sample("contemporary1")
    arg.add_sample("contemporary2")
    arg.thread_sample([0, 75], [0, 0], [4, 10])
    arg.add_sample("contemporary3")
    arg.thread_sample([0], [1], [1])
    arg.populate_children_and_roots()

    # Generates the following ARG of contemporaneously sample individuals
    #     ┊   ┃  ┊     ┊   ┃  ┊
    #     ┊   ┃  ┊10.00┊ ┏━┻┓ ┊
    #     ┊   ┃  ┊     ┊ ┃  ┃ ┊
    # 4.00┊ ┏━┻┓ ┊     ┊ ┃  ┃ ┊
    #     ┊ ┃  ┃ ┊     ┊ ┃  ┃ ┊
    # 1.00┊ ┃ ┏┻┓┊ 1.00┊ ┃ ┏┻┓┊
    # 0.00┊ 0 1 2┊ 0.00┊ 0 1 2┊
    #     0      75          100

    # Calculate volume of whole ARG and local trees
    tvC = arg_needle_lib.total_volume(arg)
    tvC1 = arg_needle_lib.local_volume(arg, 0, 25)
    tvC2 = arg_needle_lib.local_volume(arg, 25, 50)
    tvC3 = arg_needle_lib.local_volume(arg, 50, 75)
    tvC4 = arg_needle_lib.local_volume(arg, 75, 100)

    # Total volume for contemporaneous sample ARG
    assert tvC == pytest.approx(1200)

    # Local volume for genomic region [0, 25) in contemporaneous sample ARG
    assert tvC1 == pytest.approx(225)

    # Local volume for genomic region [25, 50) in contemporaneous sample ARG
    assert tvC2 == pytest.approx(225)

    # Local volume for genomic region [50, 75) in contemporaneous sample ARG
    assert tvC3 == pytest.approx(225)

    # Local volume for genomic region [75, 100) in contemporaneous sample ARG
    assert tvC4 == pytest.approx(525)

    # Add on aDNA samples
    arg.add_sample("aDNA1", 2)
    arg.thread_sample([0, 50, 75], [0, 0, 1], [5, 10, 6.5])
    arg.add_sample("aDNA2", 3)
    arg.thread_sample([0, 25], [3, 1], [10, 3.5])
    arg.add_sample("aDNA3", 7)
    arg.thread_sample([0, 25, 50, 75], [4, 0, 3, 1], [8, 10, 9, 9])
    arg.populate_children_and_roots()

    # Generates the following ARG of contemporaneously sample individuals
    # (where x denotes an aDNA sample)
    #10.00┊   ┏━━┻━━━┓ ┊   ┏━━━┻━━━┓┊    ┏━━┻━━┓ ┊ ┏━━━━┻━━━┓ ┊
    # 9.00┊   ┃      ┃ ┊   ┃       ┃┊    ┃    ┏┻┓┊ ┃      ┏━┻┓┊
    # 8.00┊   ┃     ┏┻┓┊   ┃       ┃┊    ┃    ┃ ┃┊ ┃      ┃  ┃┊
    # 7.00┊   ┃     ┃ x┊   ┃       x┊    ┃    x ┃┊ ┃      ┃  x┊
    # 6.50┊   ┃     ┃  ┊   ┃        ┊    ┃      ┃┊ ┃    ┏━┻┓  ┊
    # 6.00┊   ┃     ┃  ┊   ┃        ┊    ┃      ┃┊ ┃    ┃  ┃  ┊
    # 5.00┊ ┏━┻━┓   ┃  ┊ ┏━┻━━┓     ┊    ┃      ┃┊ ┃    ┃  ┃  ┊
    # 4.00┊ ┃ ┏━┻┓  ┃  ┊ ┃ ┏━━┻━┓   ┊ ┏━━┻━┓    ┃┊ ┃    ┃  ┃  ┊
    # 3.50┊ ┃ ┃  ┃  ┃  ┊ ┃ ┃  ┏━┻┓  ┊ ┃  ┏━┻┓   ┃┊ ┃  ┏━┻┓ ┃  ┊
    # 3.00┊ ┃ ┃  ┃  x  ┊ ┃ ┃  ┃  x  ┊ ┃  ┃  x   ┃┊ ┃  ┃  x ┃  ┊
    # 2.00┊ x ┃  ┃     ┊ x ┃  ┃     ┊ ┃  ┃      x┊ ┃  ┃    x  ┊
    # 1.00┊   ┃ ┏┻┓    ┊   ┃ ┏┻┓    ┊ ┃ ┏┻┓      ┊ ┃ ┏┻┓      ┊
    # 0.00┊ 3 0 1 2 4 5┊ 3 0 1 2 4 5┊ 0 1 2 4 5 3┊ 0 1 2 4 3 5┊
    #     0           25            50           75           100

    # Calculate volume of whole ARG and local times
    tvD = arg_needle_lib.total_volume(arg)
    tvD1 = arg_needle_lib.local_volume(arg, 0, 25)
    tvD2 = arg_needle_lib.local_volume(arg, 25, 50)
    tvD3 = arg_needle_lib.local_volume(arg, 50, 75)
    tvD4 = arg_needle_lib.local_volume(arg, 75, 100)

    # Total volume for aDNA sample ARG
    assert tvD == pytest.approx(2525)

    # Local volume for genomic region [0, 25) in aDNA sample ARG
    assert tvD1 == pytest.approx(650)

    # Local volume for genomic region [25, 50) in aDNA sample ARG
    assert tvD2 == pytest.approx(537.5)

    # Local volume for genomic region [50, 75) in aDNA sample ARG
    assert tvD3 == pytest.approx(637.5)

    # Local volume for genomic region [75, 100) in aDNA sample ARG
    assert tvD4 == pytest.approx(700)

    # Calculating aDNA additions by hand for tests that volume in two ARGs matches
    tvaDNA = 25*((3+5+1)+(5+2)+(1))+25*((3+5+1)+(0.5)+(3))+25*((7+1)+(0.5+6)+(2))+25*((4.5)+(0.5)+(2))
    tvaDNA1 = 25*((3+5+1)+(5+2)+(1))
    tvaDNA2 = 25*((3+5+1)+(0.5)+(3))
    tvaDNA3 = 25*((7+1)+(0.5+6)+(2))
    tvaDNA4 = 25*((4.5)+(0.5)+(2))

    # Difference in total volume between the two ARGs
    assert tvD - tvC == pytest.approx(1325)

    # Volume contribution from aDNA samples only
    assert tvaDNA == pytest.approx(1325)

    # Difference in local volume by genomic region for the two ARGs
    assert tvD1 - tvC1 == pytest.approx(425)
    assert tvD2 - tvC2 == pytest.approx(312.5)
    assert tvD3 - tvC3 == pytest.approx(412.5)
    assert tvD4 - tvC4 == pytest.approx(175)

    # Contribution from aDNA to local tree by genomic region
    assert tvaDNA1 == pytest.approx(425)
    assert tvaDNA2 == pytest.approx(312.5)
    assert tvaDNA3 == pytest.approx(412.5)
    assert tvaDNA4 == pytest.approx(175)
