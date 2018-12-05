import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import numpy as np


torch.manual_seed(1)


class LSTM(nn.Module):
    def __init__(self, input_dim, hidden_dim=200, alphabet_size=256):
        super(LSTM, self).__init__()
        self.input_dim = input_dim
        self.hidden_dim = hidden_dim
        self.alphabet_size = alphabet_size

        self.lstm = nn.LSTM(input_dim, hidden_dim)

        # The linear layer that maps from hidden state space to alphabet space
        self.hidden2alph = nn.Linear(hidden_dim, alphabet_size)
        self.hidden = self.init_hidden()

        self.train_initialised = False

    def init_hidden(self):
        # Before we've done anything, we dont have any hidden state.
        # The axes semantics are (num_layers, minibatch_size, hidden_dim)
        return (torch.zeros(1, 1, self.hidden_dim),
                torch.zeros(1, 1, self.hidden_dim))

    def forward(self, sequence):
        # embeds = self.word_embeddings(sentence)
        lstm_out, self.hidden = self.lstm(sequence, self.hidden)
        alphabet_space = self.hidden2alph(lstm_out.view(len(sequence), -1))
        symbol_scores = F.log_softmax(alphabet_space, dim=1)
        return symbol_scores


class LSTMOnline(object):
    def __init__(self, input_dim, hidden_dim=200, alphabet_size=256, lr=0.1):
        self.model = LSTM(input_dim, hidden_dim, alphabet_size)

        self.loss_function = nn.NLLLoss()
        self.lr = lr
        self.optimizer = optim.SGD(self.model.parameters(), lr=lr)

    def train_step(self, sequence, next_symbol):
        self.model.zero_grad()
        self.model.hidden = self.model.init_hidden()

        symbol_probs = self.model(sequence)

        loss = self.loss_function(symbol_probs, self._get_true_probs(next_symbol))
        loss.backward()
        self.optimizer.step()
        return

    def _get_true_probs(self, symbol):
        true_probs = np.zeros(shape=[self.model.input_dim], dtype=np.float32)
        true_probs[symbol] = 1.
        return true_probs

    def predict(self, sequence):
        self.model.zero_grad()
        self.model.hidden = self.model.init_hidden()

        symbol_probs = self.model(sequence)
        self.last_symbol_probs = symbol_probs
        return symbol_probs

    def train_post_predict(self, next_symbol):
        loss = self.loss_function(self.last_symbol_probs, self._get_true_probs(next_symbol))
        loss.backward()
        self.optimizer.step()
        return
