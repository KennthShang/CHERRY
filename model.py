import  torch
from    torch import nn
from    torch.nn import functional as F
from    layer import GraphConvolution
import numpy as np
from    config import args
import random

class encoder(nn.Module):
    def __init__(self, node_dim, input_dim, output_dim, num_features_nonzero):
        super(encoder, self).__init__()

        self.input_dim = input_dim
        self.output_dim = output_dim
        self.node_dim = 1474

        print('input dim:', input_dim)
        print('output dim:', output_dim)
        print('num_features_nonzero:', num_features_nonzero)

        # hidden = 64
        self.layers = nn.Sequential(GraphConvolution(self.node_dim, self.input_dim, self.input_dim, num_features_nonzero,
                                                     activation=F.relu,
                                                     dropout=0,
                                                     is_sparse_inputs=False,
                                                     learn_weight=True),

                                    #GraphConvolution(self.node_dim, self.input_dim, self.input_dim, num_features_nonzero,
                                    #                 activation=F.relu,
                                    #                 dropout=0.5,
                                    #                 is_sparse_inputs=False,
                                    #                 learn_weight=False),

                                    )


    def forward(self, inputs):
        x, support = inputs

        x = self.layers((x, support))
        
        return x[0]

    def l2_loss(self):

        layer = self.layers.children()
        layer = next(iter(layer))

        loss = None

        for p in layer.parameters():
            if loss is None:
                loss = p.pow(2).sum()
            else:
                loss += p.pow(2).sum()

        return loss


class decoder(nn.Module):
    def __init__(self, feature_dim, hidden_dim_1, hidden_dim_2):
        super(decoder, self).__init__()

        self.feature_dim = feature_dim
        self.hidden_dim_1 = hidden_dim_1
        self.hidden_dim_2 = hidden_dim_2

        self.model = nn.Sequential(
                        nn.Linear(self.feature_dim, self.hidden_dim_1),
                        nn.ReLU(),
                        nn.Linear(self.hidden_dim_1, self.hidden_dim_2),
                        nn.ReLU(),
                        nn.Linear(self.hidden_dim_2, 1)
                    )

        

    def forward(self, inputs):
        logit = self.model(inputs)
        return logit

