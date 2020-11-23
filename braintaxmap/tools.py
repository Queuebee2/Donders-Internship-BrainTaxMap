from random import random
from time import sleep

FUZZING = True


def fuzz(amt=False, fuzzing=FUZZING):
    """Inspired by Raymond Hettingers python talks
    https://github.com/rhettinger
    """
    if fuzzing:
        if not amt:
            sleep(random())
        else:
            sleep(amt)


def centr_print(title='', motif='- ', amt=60):
    """center some text with custom side motif"""
    title = title + ' ' if len(title) % 2 == 0 else title
    charsleft = abs(amt - len(title))
    side = (charsleft // 2 // len(motif))
    print(f"{(side * motif)}{motif[0]} {title} {side * motif}")


def tester(TESTS):
    """

    @param TESTS:  dict of {'program':(mode, funcname)}
    """
    if len(TESTS) == 0:
        print('No programs to run!')
        return
    centr_print(f'Testing stuff!')

    for testkey, (testmode, testfunc) in TESTS.items():
        if testmode:
            centr_print(f'testing {testkey}')
            testfunc()
            centr_print(f'Done testing {testkey}')
        else:
            print(f'testing is DISABLED for {testkey}')

    centr_print(f'All testing done')


def programrunner(PROGRAMS):
    """

    @param PROGRAMS:  dict of {'program':(mode, funcname)}
    """
    if len(PROGRAMS) == 0:
        print('No programs to run!')
        return
    centr_print(f'Running programs')

    for progkey, (mode, func) in PROGRAMS.items():
        if mode:
            centr_print(f'running {progkey}')
            func()
            centr_print(f'Done testing {progkey}')
        else:
            print(f'run is DISABLED for {progkey}')

    centr_print(f'Done running programs')
