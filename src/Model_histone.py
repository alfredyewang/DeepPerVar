from __future__ import print_function
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.optim.lr_scheduler import StepLR
import numpy as np

torch.manual_seed(0)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
np.random.seed(0)
batch_size = 256
count= 0
from torch.autograd import Variable

class AnnoSubNet(nn.Module):
    '''
    The subnetwork that is used in TFN for video and audio in the pre-fusion stage
    '''

    def __init__(self):
        super(AnnoSubNet, self).__init__()

        self.fc1 = nn.Linear(7, 2)


    def forward(self, x):

        x = torch.flatten(x, 1)
        x = self.fc1(x)
        x = F.relu(x)
        return x


class SeqSubNet(nn.Module):

    def __init__(self):
        super(SeqSubNet, self).__init__()
        self.conv1 = nn.Conv2d(4, 512, (1, 26), (1, 1))
        self.bnc1 = nn.BatchNorm2d(512)
        self.pool1 = nn.MaxPool2d((1, 13), (1, 13))
        self.dropoutc1 = nn.Dropout2d(0.25)

        self.lstm_hidden_size = 128
        self.lstm_n_layers = 4
        self.rnn = nn.LSTM(input_size=512,
                           hidden_size=self.lstm_hidden_size,
                           dropout=0.5,
                           num_layers=self.lstm_n_layers, batch_first=True, bidirectional=True)


        self.dropout1 = nn.Dropout2d(0.5)
        self.dropout2 = nn.Dropout2d(0.5)
        self.fc1 = nn.Linear(38656, 1024)
        self.bn1 = nn.BatchNorm1d(num_features=1024)
        self.fc2 = nn.Linear(1024, 128)
        self.bn2 = nn.BatchNorm1d(num_features=128)
        self.fc3 = nn.Linear(128, 64)


    def forward(self, x):
        x = self.conv1(x)
        x = self.bnc1(x)
        x = torch.relu(x)
        x = self.pool1(x)
        x = self.dropoutc1(x)

        x = x.view([-1, x.size()[1], x.size()[3]])
        x = x.permute(0, 2, 1)
        x, (h_n, h_c) = self.rnn(x)
        x = x.reshape(x.size(0), -1)
        x = self.dropout1(x)

        x = self.fc1(x)
        x = self.bn1(x)
        x = torch.relu(x)
        #
        x = self.dropout1(x)

        x = self.fc2(x)
        x = self.bn2(x)
        x = torch.relu(x)

        x = self.fc3(x)
        return x


class Net(nn.Module):
    '''
    Implements the Tensor Fusion Networks for multimodal sentiment analysis as is described in:
    Zadeh, Amir, et al. "Tensor fusion network for multimodal sentiment analysis." EMNLP 2017 Oral.
    '''

    def __init__(self):

        super(Net, self).__init__()

        self.anno_subnet = AnnoSubNet()
        self.seq_subnet = SeqSubNet()

        self.dropout1 = nn.Dropout(0.5)
        self.fc1 = nn.Linear(195, 128)
        self.bn1 = nn.BatchNorm1d(num_features=128)
        self.dropout2 = nn.Dropout(0.5)
        self.fc2 = nn.Linear(128, 64)
        self.bn2 = nn.BatchNorm1d(num_features=64)
        self.fc3 = nn.Linear(64, 1)

    def forward(self, ann_x, seq_x, batch_size):


        seq_h = self.seq_subnet(seq_x)
        anno_h = self.anno_subnet(ann_x)
        if anno_h.is_cuda:
            DTYPE = torch.cuda.FloatTensor
        else:
            DTYPE = torch.FloatTensor


        anno_h = torch.cat((Variable(torch.ones(batch_size, 1).type(DTYPE), requires_grad=False), anno_h), dim=1)
        seq_h = torch.cat((Variable(torch.ones(batch_size, 1).type(DTYPE), requires_grad=False), seq_h), dim=1)

        fusion_tensor = torch.bmm(anno_h.unsqueeze(2), seq_h.unsqueeze(1))

        fusion_tensor = torch.flatten(fusion_tensor, 1)
        x = self.dropout1(fusion_tensor)
        x = self.fc1(x)
        x = self.bn1(x)
        x = torch.relu(x)
        x = self.dropout2(x)
        x = self.fc2(x)
        x = self.bn2(x)
        x = torch.relu(x)
        x = self.fc3(x)
        x = torch.relu(x)
        return x


def test(model, device, x_test):
    model.eval()

    batch_size2 = 512
    preds = []


    x_test = x_test.astype(float).transpose((0, 2, 1))

    tmp = np.zeros((1,7))
    x_clinical = np.tile(tmp, (x_test.shape[0], 1))
    x_test = x_test.reshape(x_test.shape[0], 4, 1, 2000)
    x_test = torch.from_numpy(x_test).float()
    x_clinical = torch.from_numpy(x_clinical).float()

    with torch.no_grad():
        total_batch = int(np.shape(x_test)[0] / batch_size2)

        for i in range(total_batch + 1):
            if i == total_batch:
                x_batch, x_cli = x_test[i * batch_size2:x_test.shape[0]], \
                                        x_clinical[i * batch_size2:x_test.shape[0]]
            else:
                x_batch , x_cli= x_test[i * batch_size2:i * batch_size2 + batch_size2], \
                                        x_clinical[i * batch_size2:i * batch_size2 + batch_size2]

            data,x_cli = x_batch.to(device), x_cli.to(device)
            if data.size()[0] == 0:
                break
            output = model(x_cli,data, data.shape[0])
            preds.extend(output.data.cpu().numpy())
    preds = np.array(preds)**2
    return preds

def main(x_test, model_dir):


    torch.manual_seed(0)

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    model = Net()

    model.load_state_dict(torch.load('{}/DeepPerVar_histone.pt'.format(model_dir), map_location='cpu'))
    model = model.to(device)

    optimizer = optim.Adam(model.parameters(), lr=5e-5
                           , weight_decay=0.0001
                           )

    scheduler = StepLR(optimizer, step_size=1, gamma=0.99)

    res = test(model, device,x_test)
    return res

