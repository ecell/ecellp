import os.path
from urllib2 import Request, urlopen, URLError


def fetch_url(url, filename):
    if os.path.isfile(filename):
        # File alraedy exists
        return False

    req = Request(url)
    try:
        response = urlopen(req)
    except URLError, e:
        if hasattr(e, 'reason'):
            print 'We failed to reach a server.'
            print 'Reason: ', e.reason
            return False
        elif hasattr(e, 'code'):
            print 'The server couldn\'t fulfill the request.'
            print 'Error code: ', e.code
            return False
    else:
        # everything is fine
        root_path = os.path.dirname(filename)
        if not os.path.isdir(root_path):
            os.makedirs(root_path)
            print "Made directories [%s]" % (root_path)

        print "Dowloading [%s] from [%s] ..." % (
            filename, url)
        with open(filename, "w") as fout:
            fout.write(response.read())
    return True

