from . import tl
from . import pl
from . import dataset

# might be supressing too many future warnings.
# when in tests, delete below two lines to check future warnings.
import warnings
warnings.simplefilter('ignore', FutureWarning)

__copyright__ = 'Copyright (C) 2021 Reiichi Sugihara, Yuki Kato, Tomoya Mori, Yukio Kawahara'
__version__ = '1.0.0'
__license__ = 'BSD-3-Clause'
__author__ = 'Reiichi Sugihara, Yuki Kato, Tomoya Mori, Yukio Kawahara'
__author_email__ = 'reiichi.sugihara@gmail.com'
__url__ = 'https://github.com/ykat0/capital'
