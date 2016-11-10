import unittest

class ETCTestCase(unittest.TestCase):

    def setUp(self):
        self.browser = webdriver.Firefox()
        self.addCleanup(self.browser.quit)

    def testPageTitle(self):
        self.browser.get('http://etcbrady.stsci.edu:4242/')
        self.assertIn('DEV ETC', self.browser.title)

if __name__ == '__main__':
    unittest.main(verbosity=2)