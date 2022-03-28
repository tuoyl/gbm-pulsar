{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "124e3731-121c-4213-889f-7eba021f0356",
   "metadata": {},
   "outputs": [],
   "source": [
    "import matplotlib.pyplot as plt"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "id": "0ac2ab48-8fb1-4742-8015-8343048d6fd3",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<Figure size 432x288 with 1 Axes>"
      ]
     },
     "metadata": {
      "needs_background": "light"
     },
     "output_type": "display_data"
    }
   ],
   "source": [
    "plt.figure()\n",
    "plt.plot(1, 1)\n",
    "plt.show()"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.7.10"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
