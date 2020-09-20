import sys, argparse, json

class RPMReleases():
    def __init__(self, filename):
        self._filename = filename

    def __enter__(self):
        self.load()
        return self

    def __exit__(self, _, __, ___):
        self.dump()

    def load(self):
        try:
            with open(self._filename, 'r') as f:
                self._releases = json.load(f)
        except IOError:
            pass

    def dump(self):
        with open(self._filename, 'w+') as f:
            json.dump(self._releases, f)

    def get_next_release(self, package, version, increment = True):
        if package not in self._releases:
            self._releases[package] = {}

        if version not in self._releases[package]:
            self._releases[package][version] = 0
            return self._releases[package][version]

        if increment:
            self._releases[package][version] += 1

        return self._releases[package][version]

    def set_release(self, package, version, release):
        if package not in self._releases:
            self._releases[package] = {}

        self._releases[package][version] = release

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("jsonfile",
        help = "The json file to store versions in")
    parser.add_argument("package",
        help = "The software package for which you want a release number")
    parser.add_argument("version",
        help = "The version of the package for which you want a release number")
    parser.add_argument("--current", action = "store_true",
        help = "Don't increment any counters, just print the current value")
    parser.add_argument("--full", action = "store_true",
        help = "Print the full package name")

    args = parser.parse_args()

    with RPMReleases(args.jsonfile) as releases:
        inc = not args.current
        release = releases.get_next_release(args.package, args.version, increment = inc)

        if not args.full:
            print(release)
        else:
            print('%s-%s-%s' % (args.package, args.version, release))
    
