# Debian Packaging
To do a new release:
* Update the pkg/debian/changelog file using
  `dch -v $(<VERSION)-1 -c pkg/debian/changelog`
* When happy, finalise it using
  `dch -r -c pkg/debian/changelog`
* Commit the changelog
